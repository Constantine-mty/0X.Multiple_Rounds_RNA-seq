"""GO biological process classification rules and reusable entry points.

The classifier uses Description and optional GO ID only. GSEA-specific columns
and statistics are handled in gsea_classification.py.

Example:
    classify_go_bp_term("DNA replication", "GO:0006260")
    classify_go_bp_table(dataframe_with_ID_and_Description)
"""

import re
import pandas as pd

# ============================================================================
# GO 层级模块（可选加载）
# ============================================================================
_HIERARCHY = None
_HIERARCHY_AVAILABLE = False

def _try_load_hierarchy():
    """尝试加载 GO 层级模块，失败则退化为纯关键词模式。"""
    global _HIERARCHY, _HIERARCHY_AVAILABLE
    try:
        from .go_ontology import GOHierarchy
        _HIERARCHY = GOHierarchy()
        _HIERARCHY_AVAILABLE = True
    except Exception:
        _HIERARCHY = None
        _HIERARCHY_AVAILABLE = False

_try_load_hierarchy()

# ============================================================================
# 关键词规则
# 每条规则: (类别名, [(关键词, 匹配模式), ...])
# 匹配模式: "w" = 词边界匹配（默认），"s" = 子串匹配（用于前缀型关键词）
# 规则按优先级从高到低排列，第一条命中的规则决定类别（first match wins）
# ============================================================================

# --- 背景集标记：以下类别属于原「去分化背景项」，Background=True ---
_BACKGROUND_CATEGORIES = {
    "神经系统发育", "心脏发育", "血管与内皮发育", "骨骼与软骨发育",
    "肌肉发育", "上皮与皮肤发育", "器官发育", "胚胎与体轴发育",
    "形态发生", "细胞分化与命运决定", "系统生理与行为",
    "发育·其他（去分化背景项）", "植物发育与生理",
}

# --- 中→英类别名映射 ---
_CATEGORY_EN = {
    "线粒体功能与OXPHOS":          "Mitochondrial Function & OXPHOS",
    "DNA复制":                    "DNA Replication",
    "DNA损伤应答与修复":            "DNA Damage Response & Repair",
    "有丝分裂与染色体分离":          "Mitosis & Chromosome Segregation",
    "翻译与核糖体生物发生":          "Translation & Ribosome Biogenesis",
    "RNA加工与剪接":               "RNA Processing & Splicing",
    "转录与染色质":                 "Transcription & Chromatin",
    "脂质与固醇代谢":               "Lipid & Sterol Metabolism",
    "糖代谢与糖酵解":               "Carbohydrate Metabolism & Glycolysis",
    "氨基酸代谢":                  "Amino Acid Metabolism",
    "核苷酸与嘌呤代谢":             "Nucleotide & Purine Metabolism",
    "代谢（通用/其他）":            "General Metabolism (Other)",
    "离子转运与稳态":               "Ion Transport & Homeostasis",
    "蛋白质折叠与降解":             "Protein Folding & Degradation",
    "蛋白质定位与靶向":             "Protein Localization & Targeting",
    "蛋白质修饰与加工":             "Protein Modification & Processing",
    "蛋白质复合体组装与组织":        "Protein Complex Assembly & Organization",
    "小分子与跨膜转运":             "Small Molecule & Membrane Transport",
    "细胞器与膜组织":               "Organelle & Membrane Organization",
    "细胞壁组织与生物发生":          "Cell Wall Organization & Biogenesis",
    "光合作用与叶绿体功能":          "Photosynthesis & Chloroplast Biology",
    "植物发育与生理":               "Plant Development & Physiology",
    "微生物生长与生物膜":            "Microbial Growth & Biofilm",
    "宿主-微生物相互作用":           "Host-Microbe Interaction",
    "解毒与异生物清除":             "Detoxification & Xenobiotic Clearance",
    "病毒生命周期":                 "Viral Life Cycle",
    "刺激与激素响应（其他）":        "Stimulus & Hormone Response (Other)",
    "细胞外基质组织":               "Extracellular Matrix Organization",
    "细胞黏附":                    "Cell Adhesion",
    "细胞迁移与趋化":               "Cell Migration & Chemotaxis",
    "细胞骨架组织、连接与极性":       "Cytoskeleton Organization, Junctions & Polarity",
    "T细胞激活与信号":              "T Cell Activation & Signaling",
    "B细胞与抗体":                 "B Cell & Antibody",
    "白细胞激活与迁移":             "Leukocyte Activation & Migration",
    "细胞因子与白细胞介素响应":       "Cytokine & Interleukin Response",
    "抗原加工与递呈":               "Antigen Processing & Presentation",
    "杀伤功能":                    "Cytotoxicity & Cell Killing",
    "免疫效应与耐受":               "Immune Effector & Tolerance",
    "凋亡、程序性死亡与自噬":         "Apoptosis, Programmed Cell Death & Autophagy",
    "应激反应":                    "Stress Response",
    "干扰素与抗病毒反应":           "Interferon & Antiviral Response",
    "细胞周期调控与增殖":            "Cell Cycle Regulation & Proliferation",
    "生长因子/RTK信号":            "Growth Factor/RTK Signaling",
    "G蛋白与MAPK级联信号":         "G Protein & MAPK Cascade Signaling",
    "形态发生素信号":               "Morphogen Signaling",
    "细胞内信号转导(通用)":          "Intracellular Signal Transduction (General)",
    "囊泡运输与分泌":               "Vesicle Transport & Secretion",
    "神经系统发育":                 "Nervous System Development",
    "心脏发育":                    "Heart Development",
    "血管与内皮发育":               "Vascular & Endothelial Development",
    "骨骼与软骨发育":               "Skeletal & Cartilage Development",
    "肌肉发育":                    "Muscle Development",
    "上皮与皮肤发育":               "Epithelial & Skin Development",
    "器官发育":                    "Organ Development",
    "胚胎与体轴发育":               "Embryonic & Body Axis Development",
    "形态发生":                    "Morphogenesis",
    "细胞分化与命运决定":            "Cell Differentiation & Fate Determination",
    "伤口愈合与再生":               "Wound Healing & Regeneration",
    "系统生理与行为":               "System Physiology & Behavior",
    "发育·其他（去分化背景项）":      "Development-Other (Differentiation Background)",
    "其他/未归类":                  "Other/Uncategorized",
}

_SIGNALING_CATEGORIES = {
    "Growth Factor/RTK Signaling",
    "G Protein & MAPK Cascade Signaling",
    "Morphogen Signaling",
    "Intracellular Signal Transduction (General)",
}

_DNA_GENOME_STABILITY_CATEGORIES = {
    "DNA Replication",
    "DNA Damage Response & Repair",
    "Mitosis & Chromosome Segregation",
}

_MIGRATION_ADHESION_CYTOSKELETON_CATEGORIES = {
    "Cell Migration & Chemotaxis",
    "Cell Adhesion",
    "Cytoskeleton Organization, Junctions & Polarity",
}

_GENERAL_METABOLISM_CATEGORIES = {
    "Lipid & Sterol Metabolism",
    "Carbohydrate Metabolism & Glycolysis",
    "General Metabolism (Other)",
    "Amino Acid Metabolism",
    "Nucleotide & Purine Metabolism",
}

_IFN_ANTIVIRAL_RE = re.compile(
    r"\b(?:interferon|antiviral|anti-viral|response to (?:virus|viral)|"
    r"defense response to (?:virus|viral)|viral genome|viral replication|"
    r"rig[- ]i|mda[- ]?5)\b",
    re.IGNORECASE,
)

RULES = [
    # === 1. 能量与生物合成（高优先级）===
    ("线粒体功能与OXPHOS", [
        ("mitochondr", "s"), ("oxidative phosphorylation", "w"),
        ("electron transport", "w"), ("oxphos", "w"),
        ("atp synthesis coupled", "w"), ("respiratory chain", "w"),
        ("iron-sulfur", "w"), ("fe-s cluster", "w"),
        ("nadph dehydrogenase", "w"), ("complex iv", "w"),
        ("cytochrome c oxidase", "w"),
        ("nadh dehydrogenase", "w"), ("cytochrome complex", "w"),
        ("metallo-sulfur", "w"), ("iron-sulfur cluster", "w"),
        ("precursor metabolites", "w"), ("cellular respiration", "w"),
    ]),

    # === 2. DNA 相关（复制 > 损伤修复 > 有丝分裂）===
    ("DNA复制", [
        ("dna replication", "w"), ("replication fork", "w"),
        ("replication initiation", "w"), ("dna unwinding", "w"),
        ("dna helicase", "w"), ("dna duplex", "w"),
        ("dna conformation", "w"), ("dna geometric", "w"),
        ("dna topological", "w"), ("dna strand elongation", "w"),
        ("dna biosynthetic", "w"), ("dna synthesis", "w"),
        ("replication checkpoint", "w"),
    ]),

    ("DNA损伤应答与修复", [
        ("dna damage", "w"), ("double-strand break", "w"),
        ("replication stress", "w"), ("dna lesion", "w"),
        ("somatic hypermutation", "w"), ("dna repair", "w"),
        ("recombination", "w"), ("interstrand", "w"),
        ("fanconi", "w"), ("cross-link", "w"), ("crosslink", "w"),
        ("telomere", "w"), ("telomerase", "w"),
        ("dna strand break", "w"), ("dna crosslink", "w"),
        ("postreplication repair", "w"), ("meiotic recombination", "w"),
        ("dna modification", "w"), ("dna deamination", "w"),
        ("dna demethylation", "w"), ("dna metabolic", "w"),
        ("dna catabolic", "w"), ("dna ligation", "w"),
        ("nucleotide-excision", "w"), ("base-excision", "w"),
        ("mismatch repair", "w"),
        ("dna dealkylation", "w"), ("transposition", "w"),
        ("retrotransposition", "w"),
    ]),

    ("有丝分裂与染色体分离", [
        ("mitotic", "w"), ("mitosis", "w"), ("spindle", "w"),
        ("chromosome segregation", "w"), ("kinetochore", "w"),
        ("centromere", "w"), ("chromatid", "w"),
        ("metaphase", "w"), ("anaphase", "w"), ("telophase", "w"),
        ("prophase", "w"), ("nuclear division", "w"),
        ("cytokinesis", "w"), ("chromosome alignment", "w"),
        ("chromosome separation", "w"), ("chromosome condensation", "w"),
        ("chromosome localization", "w"), ("localization to chromosome", "w"),
        ("meiotic", "w"), ("meiosis", "w"),
        ("nuclear envelope disassembly", "w"),
        ("nuclear envelope reformation", "w"),
        ("centrosome", "w"), ("condensin", "w"), ("cohesin", "w"),
        ("sister chromatid", "w"), ("chromosome organization", "w"),
    ]),

    # === 3. 翻译与 RNA ===
    ("翻译与核糖体生物发生", [
        ("translation", "w"),
        ("ribosome", "w"), ("ribosomal", "w"),
        ("peptide biosynthesis", "w"), ("protein biosynthetic", "w"),
        ("aminoacylation", "w"),  # tRNA aminoacylation → 翻译
        ("ribosome assembly", "w"), ("ribosome biogenesis", "w"),
        ("ribosome disassembly", "w"), ("ribosome localization", "w"),
        ("ribosomal subunit", "w"),
        ("translational initiation", "w"), ("translational termination", "w"),
        ("translational elongation", "w"), ("translational fidelity", "w"),
    ]),

    ("RNA加工与剪接", [
        ("rna processing", "w"), ("splicing", "w"),
        ("spliceosome", "w"), ("spliceosomal", "w"),
        ("snrnp", "w"), ("mrna processing", "w"),
        ("mrna splicing", "w"), ("mrna metabolic", "w"),
        ("trna processing", "w"), ("trna modification", "w"),
        ("trna methylation", "w"), ("trna catabolic", "w"),
        ("trna wobble", "w"), ("trna 5'-leader", "w"),
        ("rrna processing", "w"), ("rrna modification", "w"),
        ("rrna methylation", "w"),
        ("rna modification", "w"), ("rna editing", "w"),
        ("cajal body", "w"), ("rna 3'-end", "w"),
        ("rna 5'-end", "w"), ("rna methylation", "w"),
        ("ncrna", "w"), ("rna capping", "w"),
        ("polyadenylation", "w"), ("rna localization", "w"),
        ("rna export", "w"), ("rna transport", "w"),
        ("rna splicing", "w"), ("snrna", "w"),
        ("pseudouridine", "w"), ("mrna transport", "w"),
        ("mrna export", "w"), ("mrna 3'-end", "w"),
        ("mrna modification", "w"), ("ncrna processing", "w"),
        ("ncrna metabolic", "w"), ("ncrna transcription", "w"),
        ("snrna processing", "w"), ("rna catabolic", "w"),
        ("exonucleolytic", "w"), ("deadenylated", "w"),
        ("histone mrna", "w"),
        ("prespliceosome", "w"), ("small nucleolar", "w"),
        ("mirna", "w"), ("rna surveillance", "w"),
        ("trna", "w"),
    ]),

    ("转录与染色质", [
        ("transcription", "w"), ("chromatin", "w"),
        ("nucleosome", "w"), ("histone", "w"),
        ("dna methylation", "w"), ("gene expression", "w"),
        ("rna polymerase", "w"), ("transcriptional", "w"),
        ("epigenetic", "w"), ("chromosomal organization", "w"),
        ("chromatin remodeling", "w"), ("chromatin assembly", "w"),
        ("chromatin modification", "w"), ("chromatin organization", "w"),
        ("protein-dna complex", "w"), ("retrotransposon", "w"),
    ]),

    # === 4. 代谢 ===
    ("脂质与固醇代谢", [
        ("lipid", "w"), ("fatty acid", "w"), ("cholesterol", "w"),
        ("sterol", "w"), ("phospholipid", "w"), ("steroid", "w"),
        ("sphingolipid", "w"), ("glycerolipid", "w"),
        ("eicosanoid", "w"), ("lipoprotein", "w"),
        ("glycerol", "w"), ("glycerophospholipid", "w"),
        ("triacylglycerol", "w"), ("phinganine", "w"),
    ]),

    ("糖代谢与糖酵解", [
        ("glycolysis", "w"), ("gluconeogenesis", "w"),
        ("carbohydrate", "w"), ("glucose", "w"),
        ("pyruvate", "w"), ("pentose", "w"),
        ("glycogen", "w"), ("monosaccharide", "w"),
        ("hexose", "w"), ("galactose", "w"),
        ("glycolytic", "w"), ("fructose", "w"),
        ("alditol", "w"), ("fucose", "w"),
    ]),

    ("氨基酸代谢", [
        ("amino acid", "w"), ("transamination", "w"),
        ("urea cycle", "w"), ("amino acid metabolic", "w"),
        ("glutamine", "w"), ("glutathione metabolic", "w"),
        ("glutathione", "w"), ("cysteine", "w"),
        ("methionine", "w"), ("tryptophan", "w"),
        ("tyrosine", "w"), ("phenylalanine", "w"),
        ("serine", "w"), ("glycine", "w"),
        ("branched-chain amino acid", "w"),
        ("amide catabolic", "w"), ("amino acid catabolic", "w"),
    ]),

    ("核苷酸与嘌呤代谢", [
        ("nucleotide", "w"), ("purine", "w"), ("pyrimidine", "w"),
        ("nucleobase", "w"), ("nucleoside", "w"),
        ("gmp biosynthetic", "w"), ("imp biosynthetic", "w"),
        ("amp biosynthetic", "w"), ("nucleobase biosynthetic", "w"),
        ("ctp biosynthetic", "w"), ("imp metabolic", "w"),
        ("cAMP metabolic", "w"), ("cAMP", "w"),
        ("cGMP", "w"), ("nucleotide catabolic", "w"),
        ("ribonucleoside", "w"), ("deoxyribonucleotide", "w"),
        ("deoxyribose", "w"),
    ]),

    # === 4b. 通用代谢（兜底代谢类，在特异代谢之后）===
    ("代谢（通用/其他）", [
        ("phosphate metabolic", "w"), ("phosphorus metabolic", "w"),
        ("folic acid", "w"), ("pteridine", "w"),
        ("tetrahydrofolate", "w"), ("nadh metabolic", "w"),
        ("small molecule", "w"), ("organic acid", "w"),
        ("carboxylic acid", "w"), ("ketone metabolic", "w"),
        ("ketone biosynthetic", "w"), ("aldehyde", "w"),
        ("hemoglobin", "w"), ("monocarboxylic acid", "w"),
        ("organic hydroxy", "w"), ("phenol-containing", "w"),
        ("pyridine-containing", "w"), ("olefinic", "w"),
        ("tertiary alcohol", "w"), ("secondary metabolic", "w"),
        ("nitric-oxide synthase", "w"), ("reactive nitrogen", "w"),
        ("cellular ketone", "w"), ("cellular aldehyde", "w"),
        ("alcohol biosynthetic", "w"),
    ]),

    # === 5. 离子与蛋白质稳态 ===
    ("离子转运与稳态", [
        ("ion transport", "w"),
        ("ion transmembrane", "w"), ("homeostasis", "w"),
        ("metal ion", "w"), ("vitamin", "w"),
        ("cofactor", "w"), ("ion homeostasis", "w"),
        ("ph homeostasis", "w"), ("osmotic", "w"),
        ("water transport", "w"), ("ion channel", "w"),
        ("metal ion homeostasis", "w"), ("potassium ion", "w"),
        ("sodium ion", "w"), ("calcium ion", "w"),
        ("iron ion", "w"), ("zinc ion", "w"),
        ("copper ion", "w"), ("magnesium", "w"),
        ("chloride transport", "w"), ("solute transport", "w"),
        ("import across", "w"), ("membrane potential", "w"),
        ("voltage-gated", "w"),
        ("cation channel", "w"), ("transporter activity", "w"),
    ]),

    ("蛋白质折叠与降解", [
        ("protein folding", "w"), ("protein unfolding", "w"), ("proteolysis", "w"),
        ("ubiquitin", "w"), ("proteasome", "w"),
        ("protein catabolic", "w"),
        ("chaperone", "w"),
        ("protein modification", "w"), ("peptidase", "w"),
        ("endopeptidase", "w"),
        ("protease", "w"), ("deubiquitination", "w"),
        ("protein refolding", "w"), ("protein maturation", "w"),
        ("protein glycosylation", "w"), ("protein lipidation", "w"),
        ("protein alkylation", "w"), ("protein acylation", "w"),
        ("protein autoprocessing", "w"), ("protein methylation", "w"),
        ("protein hydroxylation", "w"), ("protein phosphorylation", "w"),
        ("protein ubiquitination", "w"), ("protein dephosphorylation", "w"),
        ("protein acetylation", "w"), ("protein sumoylation", "w"),
        ("protein neddylation", "w"), ("protein glutathionylation", "w"),
        ("protein polymerization", "w"), ("protease inhibitor", "w"),
        ("protein processing", "w"),
        ("zymogen", "w"), ("peptide cross-linking", "w"),
        ("deacetylation", "w"),
        ("k11-linked", "w"), ("protein deacetylation", "w"),
        ("glycoprotein", "w"),
    ]),

    # === 6. 细胞外基质（新增）===
    ("细胞外基质组织", [
        ("extracellular matrix", "w"), ("extracellular structure", "w"),
        ("collagen", "w"), ("fibril organization", "w"),
        ("basement membrane", "w"), ("ecm", "w"),
        ("external encapsulating structure", "w"),
        ("glycosaminoglycan", "w"), ("proteoglycan", "w"),
        ("mucopolysaccharide", "w"), ("aminoglycan", "w"),
    ]),

    # === 7. 细胞结构/运动 ===
    ("细胞黏附", [
        ("cell adhesion", "w"), ("adhesion", "w"),
        ("cadherin", "w"), ("integrin", "w"),
        ("cell-cell", "w"), ("homophilic", "w"),
        ("cell-substrate", "w"), ("focal adhesion", "w"),
        ("cell-matrix", "w"), ("cell adhesion molecule", "w"),
        ("leukocyte tethering", "w"), ("cell recognition", "w"),
    ]),

    ("细胞迁移与趋化", [
        ("migration", "w"), ("chemotaxis", "w"),
        ("locomotion", "w"), ("taxis", "w"),
        ("cell migration", "w"), ("motility", "w"),
    ]),

    ("细胞骨架组织、连接与极性", [
        ("cytoskeleton", "w"), ("actin", "w"),
        ("microtubule", "w"), ("filament", "w"),
        ("cytoskeleton organization", "w"), ("junction", "w"),
        ("polarity", "w"), ("cell projection", "w"),
        ("cilium", "w"), ("cilia", "w"), ("flagellum", "w"),
        ("lamellipodium", "w"), ("filopodium", "w"),
        ("ruffle", "w"), ("leading edge", "w"),
        ("cell shape", "w"), ("microtubule organizing", "w"),
        ("cellular protrusion", "w"), ("actomyosin", "w"),
        ("supramolecular fiber", "w"), ("fiber organization", "w"),
        ("centriole", "w"), ("bleb", "w"),
    ]),

    # === 8. 免疫（T/B/白细胞 在 通用增殖/背景 之前）===
    ("T细胞激活与信号", [
        ("t cell", "w"), ("alpha-beta t cell", "w"),
        ("t cell activation", "w"), ("t cell receptor", "w"),
        ("t cell proliferation", "w"), ("cd4-positive", "w"),
        ("cd8-positive", "w"), ("thymocyte", "w"),
        ("t cell selection", "w"), ("alpha-beta t", "w"),
        ("gamma-delta t", "w"), ("t-helper", "w"),
        ("t cell differentiation", "w"), ("t cell tolerance", "w"),
        ("t cell apoptotic", "w"), ("t cell homeostasis", "w"),
        ("regulatory t cell", "w"), ("activated t cell", "w"),
        ("t cell lineage", "w"),
    ]),

    ("B细胞与抗体", [
        ("b cell", "w"), ("immunoglobulin", "w"),
        ("isotype", "w"), ("b cell receptor", "w"),
        ("b cell proliferation", "w"), ("b cell activation", "w"),
        ("antibody", "w"), ("b cell differentiation", "w"),
    ]),

    ("白细胞激活与迁移", [
        ("leukocyte activation", "w"), ("leukocyte migration", "w"),
        ("leukocyte chemotaxis", "w"), ("myeloid", "w"),
        ("neutrophil", "w"), ("macrophage", "w"),
        ("phagocytosis", "w"), ("mast cell", "w"),
        ("granulocyte", "w"), ("monocyte", "w"),
        ("dendritic", "w"), ("leukocyte proliferation", "w"),
        ("myeloid leukocyte", "w"), ("astrocyte activation", "w"),
        ("microglial", "w"),
        ("extravasation", "w"),
    ]),

    ("细胞因子与白细胞介素响应", [
        ("cytokine", "w"), ("interleukin", "w"),
        ("cytokine production", "w"), ("cytokine response", "w"),
        ("interleukin production", "w"), ("lymphokine", "w"),
        ("chemokine", "w"), ("response to chemokine", "w"),
        ("cellular response to chemokine", "w"),
        ("tnf", "w"), ("tumor necrosis factor", "w"),
        ("colony-stimulating factor", "w"),
        ("growth factor stimulus", "w"),
    ]),

    # === 9. 杀伤与免疫效应（杀伤 在 免疫效应 之前，两者在背景之前）===
    ("杀伤功能", [
        ("cell killing", "w"), ("leukocyte mediated immunity", "w"),
        ("cytotoxicity", "w"), ("regulation of killing", "w"),
        ("target cell", "w"), ("leukocyte mediated cytotoxicity", "w"),
        ("killing of cells", "w"),
    ]),

    ("免疫效应与耐受", [
        ("immune response", "w"), ("immunity", "w"),
        ("tolerance", "w"), ("tolerance induction", "w"),
        ("nk cell", "w"), ("natural killer", "w"),
        ("humoral", "w"), ("complement", "w"),
        ("defense response to bacterium", "w"), ("antimicrobial", "w"),
        ("antibacterial", "w"), ("inflammatory", "w"),
        ("immunological", "w"), ("immune system process", "w"),
        ("immune ", "w"), ("hemopoiesis", "w"),
        ("hematopoietic", "w"), ("defense response to protozoan", "w"),
        ("defense response", "w"), ("immune receptor", "w"),
        ("antigen receptor", "w"), ("lymphocyte", "w"),
        ("leukocyte", "w"), ("myeloid cell", "w"),
        ("symbiont", "w"), ("host", "w"),
        ("interspecies", "w"), ("antibiotic", "w"),
        ("biotic stimulus", "w"), ("response to bacterium", "w"),
        ("molecule of bacterial origin", "w"),
        ("symbiotic", "w"), ("intraspecies", "w"),
        ("modulation of process", "w"), ("inflammasome", "w"),
        ("regulation of cell activation", "w"),
    ]),

    # === 10. 死亡与应激 ===
    ("凋亡、程序性死亡与自噬", [
        ("apoptosis", "w"), ("programmed cell death", "w"),
        ("cell death", "w"), ("necroptosis", "w"),
        ("pyroptosis", "w"), ("anoikis", "w"),
        ("apoptotic", "w"), ("autophagy", "w"),
        ("macroautophagy", "w"), ("autophagosome", "w"),
        ("chaperone-mediated", "w"), ("apoptotic signaling", "w"),
        ("regulation of cell death", "w"),
    ]),

    ("干扰素与抗病毒反应", [
        ("interferon", "w"), ("antiviral", "w"),
        ("anti-viral", "w"),
        ("rig-i", "w"), ("rig i", "w"),
        ("mda-5", "w"), ("mda5", "w"),
        ("response to virus", "w"), ("response to viral", "w"),
        ("viral genome", "w"), ("viral replication", "w"),
        ("type i interferon", "w"),
        ("type ii interferon", "w"), ("defense response to virus", "w"),
    ]),

    ("应激反应", [
        ("stress", "w"), ("heat shock", "w"),
        ("unfolded protein", "w"), ("hypoxia", "w"),
        ("oxidative stress", "w"), ("response to uv", "w"),
        ("response to radiation", "w"), ("response to drug", "w"),
        ("response to xenobiotic", "w"), ("response to oxidative", "w"),
        ("response to heat", "w"), ("response to cold", "w"),
        ("response to osmotic", "w"), ("response to starvation", "w"),
        ("response to nutrient", "w"), ("response to metal", "w"),
        ("response to oxygen", "w"), ("response to light", "w"),
        ("response to ph", "w"), ("response to lipopolysaccharide", "w"),
        ("response to injury", "w"), ("response to acid", "w"),
        ("response to peptide", "w"), ("response to wounding", "w"),
        ("genotoxic", "w"), ("reactive oxygen", "w"),
        ("nitric oxide", "w"), ("antioxidant", "w"),
        ("xenobiotic", "w"), ("circadian", "w"),
        ("response to ethanol", "w"), ("response to ketone", "w"),
        ("response to glucocorticoid", "w"), ("response to retinoic acid", "w"),
        ("response to dopamine", "w"), ("response to activity", "w"),
        ("response to cAMP", "w"), ("response to calcium", "w"),
        ("response to monoamine", "w"), ("response to catecholamine", "w"),
        ("response to gonadotropin", "w"),
        ("response to organic cyclic", "w"),
        ("response to organonitrogen", "w"),
        ("response to nitrogen compound", "w"),
        ("hydrogen peroxide", "w"),
        ("response to prostaglandin", "w"),
        ("response to corticosteroid", "w"),
        ("response to organophosphorus", "w"),
        ("response to external stimulus", "w"),
        ("response to amphetamine", "w"),
        ("response to decreased oxygen", "w"),
        ("response to amine", "w"),
        ("response to mechanical", "w"),
        ("response to amyloid", "w"),
        ("response to inorganic", "w"),
        ("response to estrogen", "w"),
        ("response to insulin", "w"),
        ("response to manganese", "w"),
        ("response to alcohol", "w"),
        ("response to dietary", "w"),
        ("response to estradiol", "w"),
        ("response to angiotensin", "w"),
        ("response to extracellular", "w"),
        ("response to lectin", "w"),
        ("response to protozoan", "w"),
        ("response to tumor cell", "w"),
        ("heat generation", "w"),
        ("amyloid-beta", "w"),
        ("decreased oxygen", "w"),
    ]),

    # === 11. 增殖（通用，在免疫亚类之后）===
    ("细胞周期调控与增殖", [
        ("cell cycle", "w"), ("checkpoint", "w"),
        ("g1/s", "w"), ("g2/m", "w"),
        ("cell cycle transition", "w"), ("proliferation", "w"),
        ("centrosome cycle", "w"), ("stem cell proliferation", "w"),
        ("cell proliferation", "w"), ("cell cycle phase", "w"),
        ("cell growth", "w"), ("cell division", "w"),
        ("cell size", "w"), ("cellular component size", "w"),
        ("regulation of cell size", "w"),
        ("regulation of growth", "w"), ("maintenance of cell number", "w"),
    ]),

    # === 12. 信号转导 ===
    ("生长因子/RTK信号", [
        ("growth factor", "w"), ("erbb", "w"),
        ("receptor tyrosine", "w"), ("vegf", "w"),
        ("fgf", "w"), ("egf", "w"), ("pdgf", "w"),
        ("tyrosine kinase receptor", "w"),
        ("growth factor receptor", "w"),
        ("insulin receptor", "w"), ("insulin-like", "w"),
        ("leukemia inhibitory factor", "w"),
    ]),

    ("G蛋白与MAPK级联信号", [
        ("gpcr", "w"), ("g-protein", "w"), ("g protein", "w"),
        ("gtpase", "w"), ("second messenger", "w"),
        ("mapk", "w"), ("erk", "w"), ("jnk", "w"),
        ("p38", "w"), ("kinase cascade", "w"),
        ("map kinase", "w"), ("phosphatidylinositol", "w"),
        ("ras", "w"), ("rho", "w"), ("rac", "w"),
        ("cdc42", "w"), ("gtpase activity", "w"),
        ("stat", "w"), ("jak-stat", "w"),
        ("nuclear receptor", "w"), ("intracellular receptor", "w"),
        ("ppar", "w"), ("peroxisome proliferator", "w"),
        ("nmda receptor", "w"), ("neurotransmitter receptor", "w"),
    ]),

    ("形态发生素信号", [
        ("notch", "w"), ("wnt", "w"), ("bmp", "w"),
        ("hedgehog", "w"), ("smoothened", "w"),
        ("eph", "w"), ("ephrin", "w"),
        ("tgf", "w"), ("beta-catenin", "w"),
        ("wnt signaling", "w"), ("bmp signaling", "w"),
        ("tgf-beta", "w"),
    ]),

    ("细胞内信号转导(通用)", [
        ("signal transduction", "w"), ("signaling", "w"),
        ("signal", "w"), ("kinase", "w"),
        ("phosphorylation", "w"), ("cascade", "w"),
        ("hormone signaling", "w"), ("receptor signaling", "w"),
        ("transduction", "w"), ("cell-cell signaling", "w"),
        ("protein kinase", "w"),
        ("transferase activity", "w"), ("lyase activity", "w"),
        ("monooxygenase activity", "w"), ("nuclease activity", "w"),
        ("dna binding", "w"), ("regulation of binding", "w"),
        ("cell communication", "w"),
    ]),

    # === 13. 运输与定位 ===
    ("囊泡运输与分泌", [
        ("vesicle", "w"), ("exocytosis", "w"),
        ("secretion", "w"), ("endosome", "w"),
        ("golgi", "w"), ("endocytosis", "w"),
        ("lysosome", "w"), ("endomembrane", "w"),
        ("vacuole", "w"),
        ("transport vesicle", "w"),
        ("membrane invagination", "w"),
        ("membrane raft", "w"),
        ("exocytic", "w"), ("secretory granule", "w"),
        ("cytoplasmic granule", "w"),
        ("pinocytosis", "w"),
    ]),

    # === 14. 背景亚类（发育/神经/器官发生细分，近末尾但在兜底之前）===
    ("神经系统发育", [
        ("neuron", "w"), ("axon", "w"), ("synap", "s"),
        ("neural", "w"), ("glia", "w"), ("gliogenesis", "w"),
        ("dendrite", "w"), ("neurogen", "s"),
        ("brain", "w"), ("nervous", "w"),
        ("neuron projection", "w"), ("neuromuscular", "w"),
        ("neuron differentiation", "w"), ("schwann", "w"),
        ("ganglion", "w"), ("cerebellar", "w"), ("cerebellum", "w"),
        ("cerebral cortex", "w"), ("cranial nerve", "w"),
        ("facial nerve", "w"), ("action potential", "w"),
        ("neurotransmitter", "w"), ("synaptic transmission", "w"),
        ("postsynapse", "w"), ("presynapse", "w"),
        ("neural crest", "w"), ("axon guidance", "w"),
        ("axonogenesis", "w"), ("axon development", "w"),
        ("axon extension", "w"), ("axon regeneration", "w"),
        ("synapse organization", "w"), ("myelination", "w"),
        ("nervous system", "w"),
        ("fasciculation", "w"), ("ensheathment", "w"),
        ("cognition", "w"), ("sensory perception", "w"),
        ("membrane depolarization", "w"), ("membrane repolarization", "w"),
        ("collateral sprouting", "w"), ("glial cell activation", "w"),
        ("neuroinflammatory", "w"), ("myelin assembly", "w"),
    ]),

    ("心脏发育", [
        ("heart", "w"), ("cardiac", "w"),
        ("ventric", "w"), ("myocard", "w"),
        ("atrial", "w"), ("endocard", "w"),
        ("aortic valve", "w"), ("heart morphogenesis", "w"),
    ]),

    ("血管与内皮发育", [
        ("angiogenesis", "w"), ("vascul", "s"),
        ("endothel", "s"), ("blood vessel", "w"),
        ("vasculogenesis", "w"), ("aorta", "w"),
        ("artery", "w"), ("sprouting angiogenesis", "w"),
        ("blood circulation", "w"), ("coagulation", "w"),
        ("hemostasis", "w"), ("blood pressure", "w"),
        ("vasoconstriction", "w"), ("vasodilation", "w"),
        ("circulatory", "w"),
    ]),

    ("骨骼与软骨发育", [
        ("bone", "w"), ("cartilage", "w"),
        ("skeletal", "w"), ("ossif", "s"),
        ("osteob", "s"), ("chondrocyte", "w"),
        ("biomineral", "w"), ("skeletal system", "w"),
    ]),

    ("肌肉发育", [
        ("muscle", "w"), ("sarcomere", "w"),
        ("myoblast", "w"), ("myofibril", "w"),
        ("myocardial", "w"), ("striated muscle", "w"),
        ("muscle cell", "w"), ("muscle organ", "w"),
        ("muscle contraction", "w"), ("cardiac muscle", "w"),
    ]),

    ("上皮与皮肤发育", [
        ("epithel", "s"), ("epiderm", "s"),
        ("skin", "w"), ("keratin", "w"),
        ("epithelial cell", "w"), ("epithelial tube", "w"),
        ("epithelial to mesenchymal", "w"),
        ("establishment of skin barrier", "w"),
        ("hair follicle", "w"), ("hair cycle", "w"),
        ("keratinization", "w"),
    ]),

    ("器官发育", [
        ("kidney", "w"), ("nephron", "w"), ("renal", "w"),
        ("eye", "w"), ("ear", "w"), ("cochlea", "w"),
        ("lung", "w"), ("alveol", "s"),
        ("liver", "w"), ("hepat", "s"),
        ("pancrea", "s"), ("salivary", "w"),
        ("gland", "w"), ("tooth", "w"), ("palate", "w"),
        ("digestive", "w"), ("respiratory", "w"),
        ("urogenital", "w"), ("gonad", "w"),
        ("sperm", "w"), ("oocyte", "w"),
        ("fertilization", "w"), ("embryo implantation", "w"),
        ("endocrine", "w"), ("exocrine", "w"),
        ("hepaticobiliary", "w"), ("labyrinthine", "w"),
        ("glomerulus", "w"), ("endocardial cushion", "w"),
        ("appendage", "w"), ("limb", "w"),
        ("face", "w"), ("facial", "w"), ("cranial", "w"),
        ("sensory organ", "w"), ("sensory system", "w"),
        ("visual system", "w"), ("lens fiber", "w"),
        ("organ formation", "w"), ("animal organ", "w"),
        ("connective tissue", "w"), ("adipose tissue", "w"),
        ("genitalia", "w"), ("female gonad", "w"),
        ("sertoli", "w"),
        ("organ growth", "w"), ("organ induction", "w"),
        ("enamel", "w"), ("oogenesis", "w"),
    ]),

    ("胚胎与体轴发育", [
        ("embryonic", "w"), ("gastrulation", "w"),
        ("germ layer", "w"), ("blastocyst", "w"),
        ("endoderm", "w"), ("mesoderm", "w"),
        ("ectoderm", "w"), ("epiboly", "w"),
        ("axis specification", "w"), ("pattern specification", "w"),
        ("pattern formation", "w"), ("segmentation", "w"),
        ("dorsal", "w"), ("anterior/posterior", "w"),
        ("left/right", "w"), ("bilateral symmetry", "w"),
        ("blastoderm", "w"), ("embryo development", "w"),
        ("germ ring", "w"), ("notochord", "w"),
        ("developmental induction", "w"),
        ("formation of primary germ layer", "w"),
        ("endoderm formation", "w"), ("endodermal cell", "w"),
        ("axis elongation", "w"), ("tripartite", "w"),
    ]),

    ("形态发生", [
        ("morphogenesis", "w"), ("anatomical structure", "w"),
        ("branching", "w"), ("invagination", "w"),
        ("tube morphogenesis", "w"), ("tube development", "w"),
        ("epithelial tube", "w"),
        ("tube size", "w"), ("tube diameter", "w"),
        ("tube closure", "w"),
    ]),

    ("细胞分化与命运决定", [
        ("differentiation", "w"), ("cell fate", "w"),
        ("lineage commitment", "w"), ("lineage", "w"),
        ("determination", "w"), ("specification", "w"),
        ("fate commitment", "w"), ("fate determination", "w"),
        ("cell fate specification", "w"), ("cell fate determination", "w"),
        ("stem cell", "w"), ("cell maturation", "w"),
        ("developmental maturation", "w"),
        ("anatomical structure maturation", "w"),
        ("cell development", "w"),
    ]),

    ("伤口愈合与再生", [
        ("wound", "w"), ("regeneration", "w"),
        ("wound healing", "w"),
        ("tissue remodeling", "w"),
    ]),

    ("系统生理与行为", [
        ("renal", "w"), ("excretion", "w"),
        ("renal filtration", "w"), ("renal absorption", "w"),
        ("renal system", "w"), ("glomerular filtration", "w"),
        ("body fluid", "w"), ("blood pressure", "w"),
        ("behavior", "w"), ("locomotory behavior", "w"),
        ("learning", "w"), ("memory", "w"),
        ("reproduction", "w"), ("reproduct", "s"),
        ("pregnancy", "w"), ("female pregnancy", "w"),
        ("ovulation", "w"), ("gamet", "s"),
        ("contraction", "w"), ("circadian", "w"),
        ("metamorphosis", "w"), ("larval", "w"),
        ("postembryonic", "w"), ("molting", "w"),
        ("aging", "w"), ("regulation of urine", "w"),
        ("symbiont", "w"), ("host", "w"),
        ("disruption of", "w"), ("killing of cells of another", "w"),
        ("locomotory", "w"), ("balance", "w"),
        ("pigment", "w"), ("pigmentation", "w"),
        ("cellular pigment", "w"),
        ("rhythmic process", "w"), ("fibrinolysis", "w"),
        ("multi-multicellular", "w"),
    ]),

    # === 15. 发育兜底（不匹配任何背景亚类的发育术语）===
    ("发育·其他（去分化背景项）", [
        ("development", "w"), ("developmental", "w"),
        ("regionalization", "w"),
        ("boundary", "w"), ("dorsal/ventral", "w"),
        ("proximal/distal", "w"),
    ]),

    # === 16. 兜底 ===
    ("其他/未归类", []),
]

# ============================================================================
# 5 主题追踪规则（独立于类别，更细粒度）
# ============================================================================
THEME_RULES = [
    ("DNA Repair", ["dna repair", "double-strand break repair", "single strand break repair",
                    "nucleotide-excision repair", "base-excision repair", "mismatch repair",
                    "postreplication repair", "interstrand cross-link repair"]),
    ("DNA Damage", ["dna damage", "double-strand break", "dna damage response",
                    "dna damage checkpoint", "replication stress",
                    "dna strand break", "somatic hypermutation", "dna lesion"]),
    ("Cell Killing", ["cell killing", "leukocyte mediated immunity", "regulation of killing",
                      "cytotoxicity", "target cell", "killing of cells"]),
    ("Proliferation", ["proliferation", "cell cycle", "mitotic", "cell division",
                       "chromosome segregation", "spindle", "cytokinesis",
                       "centromere", "kinetochore", "condensin", "cohesin",
                       "dna replication", "mitosis", "nuclear division",
                       "g1/s", "g2/m", "chromatid", "metaphase", "anaphase",
                       "telophase", "prophase", "centrosome"]),
    ("Programmed Cell Death", ["apoptosis", "programmed cell death", "cell death",
                               "necroptosis", "pyroptosis", "anoikis", "apoptotic"]),
]

# ============================================================================
# 预编译正则
# ============================================================================
def _compile_rules(rules):
    """将规则编译为 [(category, [(compiled_regex, kw_text), ...])]"""
    compiled = []
    for cat, kws in rules:
        items = []
        for kw_spec in kws:
            if isinstance(kw_spec, tuple):
                kw, mode = kw_spec
            else:
                kw, mode = kw_spec, "w"
            if mode == "s":
                # 子串匹配
                pat = re.compile(re.escape(kw), re.IGNORECASE)
            else:
                # 词边界匹配
                pat = re.compile(r'\b' + re.escape(kw) + r'\b', re.IGNORECASE)
            items.append((pat, kw))
        compiled.append((cat, items))
    return compiled

_RULES_COMPILED = _compile_rules(RULES)

# Specific processes precede broad first-match keywords and hierarchy anchors.
_PRIORITY_RULES = [
    ("抗原加工与递呈", re.compile(r"\b(?:antigen (?:processing|presentation)|presentation of (?:exogenous|endogenous|peptide) antigen|mhc (?:class [i]+ )?(?:antigen )?presentation)\b", re.I)),
    ("伤口愈合与再生", re.compile(r"\b(?:wound(?:ing)?|tissue (?:repair|regeneration)|skin regeneration)\b", re.I)),
    ("DNA损伤应答与修复", re.compile(r"\b(?:single[ -]strand break repair|double[ -]strand break repair|removal of nonhomologous ends)\b", re.I)),
    ("应激反应", re.compile(r"\b(?:response to [\w -]*stress|stress response|unfolded protein response)\b", re.I)),
    ("光合作用与叶绿体功能", re.compile(r"\b(?:photosynth\w*|photosystem\w*|chloroplast\w*|chlorophyll\w*|thylakoid\w*|plastid\w*|photophosphorylation)\b", re.I)),
    ("植物发育与生理", re.compile(r"\b(?:floral|flower(?:ing)?|meristem|pollen|trichome|phloem|xylem|casparian strip|haustorium|seed germination|fruit (?:development|ripening)|root (?:development|growth)|leaf (?:development|senescence)|shoot (?:development|growth)|stomatal)\b", re.I)),
    ("微生物生长与生物膜", re.compile(r"\b(?:biofilm|sporulation|spore germination|spore wall|aplanospore|ascospore|oospore|prospore|akinete|germ tube|pseudohyphal|hyphal growth|hypha formation|pilus (?:assembly|retraction|biogenesis))\b", re.I)),
    ("细胞壁组织与生物发生", re.compile(r"\b(?:cell wall|peptidoglycan|lipopolysaccharide biosynth\w*)\b", re.I)),
    ("病毒生命周期", re.compile(r"\b(?:viral|virus (?:replication|entry|assembly|budding|release|life cycle)|virion|capsid)\b", re.I)),
    ("宿主-微生物相互作用", re.compile(r"\b(?:symbiont-mediated|symbiont-induced|symbiont entry|host-seeking|entry of bacterium into host|exit from host cell)\b", re.I)),
    ("解毒与异生物清除", re.compile(r"\b(?:detoxification|xenobiotic (?:catabolic|metabolic|clearance)|drug clearance)\b", re.I)),
    ("蛋白质定位与靶向", re.compile(r"\b(?:protein (?:localization|targeting|import|export|translocation|insertion|transport|transmembrane transport)|(?:targeting|localization|transport|import|export) of protein|protein-containing complex localization)\b", re.I)),
    ("蛋白质修饰与加工", re.compile(r"\b(?:protein (?:modification|maturation|processing|glycosylation|lipidation|acylation|methylation|hydroxylation|phosphorylation|autophosphorylation|dephosphorylation|acetylation|deacetylation|sumoylation|neddylation|glutathionylation|palmitoylation|prenylation|polyglycylation|deglutarylation|oligomerization)|peptidyl-[\w-]+|post-translational protein modification)\b", re.I)),
    ("RNA加工与剪接", re.compile(r"\b(?:rna repair|(?:mrna|rrna|trna|sirna|mirna) (?:stabili\w*|decapping|splice|processing|methylation)|(?:mrna|rrna) poly\(a\) tail)\b", re.I)),
    ("蛋白质复合体组装与组织", re.compile(r"\bprotein-containing complex (?:organization|assembly|disassembly|remodeling)\b", re.I)),
    ("细胞器与膜组织", re.compile(r"\b(?:(?:organelle|membrane|nucleus|nucleolus|nuclear pore|peroxisome|endoplasmic reticulum|cellular component|nuclear envelope) (?:organization|assembly|biogenesis|maintenance|disassembly|reassembly|fusion|fission)|endoplasmic reticulum (?:tubular|cisternal) network (?:organization|maintenance|formation)|endoplasmic reticulum polarization|membrane repair)\b", re.I)),
    ("小分子与跨膜转运", re.compile(r"\b(?:small molecule transport|transmembrane transport|transport of (?:small molecule|metabolite|organic|carbohydrate|lipid)|(?:drug|xenobiotic|methotrexate|glucose|sugar|amino acid|nucleotide|lipid|fatty acid|phosphate|organic acid) transport)\b", re.I)),
    ("心脏发育", re.compile(r"\bendocardium\b", re.I)),
    ("心脏发育", re.compile(r"\bcardiac muscle cell action potential\b", re.I)),
    ("肌肉发育", re.compile(r"\bskeletal muscle\b", re.I)),
    ("细胞迁移与趋化", re.compile(r"\bchemokinesis\b", re.I)),
    ("神经系统发育", re.compile(r"\b(?:nerve (?:formation|maturation|development)|photoreceptor (?:organization|development)|visual perception|sensory processing|synaptic receptor clustering)\b", re.I)),
    ("胚胎与体轴发育", re.compile(r"\bspemann organizer\b", re.I)),
    ("B细胞与抗体", re.compile(r"\bgerminal center (?:formation|organization)\b", re.I)),
    ("血管与内皮发育", re.compile(r"\bplatelet (?:activation|aggregation|degranulation)\b", re.I)),
    ("转录与染色质", re.compile(r"\b(?:heterochromatin formation|transposable element silencing)\b", re.I)),
    ("系统生理与行为", re.compile(r"\bdetection of (?:external|abiotic|light|sound|odor|taste|temperature|oxygen|carbon dioxide|humidity|dietary)\b", re.I)),
    ("系统生理与行为", re.compile(r"\b(?:penile erection|vasomotion)\b", re.I)),
    ("G蛋白与MAPK级联信号", re.compile(r"\b(?:guanylate|adenylate) cyclase activity\b", re.I)),
    ("杀伤功能", re.compile(r"\bgranzyme\b", re.I)),
    ("免疫效应与耐受", re.compile(r"\b(?:hypersensitivity|opsonization|detection of (?:bacteri\w*|virus|viral|fung\w*|yeast|protozoan|oomycete|pathogen|symbiont|tumor cell))\b", re.I)),
    ("免疫效应与耐受", re.compile(r"\btoll-like receptor\b", re.I)),
    ("囊泡运输与分泌", re.compile(r"\b(?:receptor (?:internalization|recycling)|ubiquitin-dependent endocytosis)\b", re.I)),
]

# Name-only patterns cannot reliably distinguish these specific GO processes.
_GO_ID_OVERRIDES = {
    "GO:0001778": "细胞器与膜组织",  # plasma membrane repair
    "GO:0030091": "蛋白质折叠与降解",  # protein repair
    "GO:0010206": "光合作用与叶绿体功能",  # photosystem II repair
    "GO:0110052": "代谢（通用/其他）",  # toxic metabolite repair
    "GO:0061057": "宿主-微生物相互作用",  # peptidoglycan recognition protein signaling
    "GO:0039640": "病毒生命周期",  # viral release via host peptidoglycan suppression
    "GO:0006335": "转录与染色质",  # replication-dependent chromatin assembly
}

_METABOLISM_FALLBACK = re.compile(r"\b(?:metaboli\w*|cataboli\w*|biosynth\w*|fermentation)\b", re.I)
_TRANSPORT_FALLBACK = re.compile(r"\b(?:transport|import|export|uptake|efflux|localization)\b", re.I)
_ORGANELLE_FALLBACK = re.compile(r"\b(?:(?:organelle|membrane|nucleus|nucleolus|nuclear pore|peroxisome|endoplasmic reticulum|cellular component) (?:organization|assembly|biogenesis|maintenance|disassembly)|membrane repair)\b", re.I)
_RESPONSE_FALLBACK = re.compile(r"\b(?:cellular )?response to\b", re.I)

# ============================================================================
# 核心函数
# ============================================================================

def _keyword_match(desc):
    """纯关键词匹配（词边界），返回类别名或 None。"""
    d = desc.lower()
    for cat, items in _RULES_COMPILED:
        if not items:
            return cat
        for pat, kw in items:
            if pat.search(d):
                return cat
    return None


def _hierarchy_match(go_id):
    """GO 层级锚点匹配，返回 (类别, 分数) 或 None。"""
    if not _HIERARCHY_AVAILABLE or not go_id:
        return None
    scores = _HIERARCHY.match_anchors(go_id)
    if not scores:
        return None
    best_cat = max(scores, key=scores.get)
    return (best_cat, scores[best_cat])


def _category_en(cat_cn):
    """将内部中文类别名转换为英文类别名。"""
    return _CATEGORY_EN.get(cat_cn, "Other/Uncategorized")


def assign_fine_category_cn(desc, go_id=None):
    """
    将一条 GO 描述归入中文功能类别（混合分类：关键词 + GO层级）。

    Parameters
    ----------
    desc : str
        GO 通路描述（英文）
    go_id : str, optional
        GO ID（如 "GO:0043039"），用于层级辅助分类

    Returns
    -------
    str : 中文类别名
    """
    desc = "" if pd.isna(desc) else re.sub(r"^\s*obsolete\s+", "", str(desc), flags=re.I)
    if go_id in _GO_ID_OVERRIDES:
        return _GO_ID_OVERRIDES[go_id]

    for category, pattern in _PRIORITY_RULES:
        if pattern.search(desc):
            if category == "应激反应" and _IFN_ANTIVIRAL_RE.search(desc):
                return "干扰素与抗病毒反应"
            if category == "病毒生命周期" and re.search(r"\b(?:response to (?:virus|viral)|defense response|antiviral|interferon)\b", desc, re.I):
                continue
            if category == "小分子与跨膜转运" and re.search(r"\b(?:ion|cation|anion|sodium|potassium|calcium|chloride|magnesium|zinc|copper|iron|proton|fluoride|iodide|sulfite|sulfate|bicarbonate|ammonium|nitrate|nitrite)\b", desc, re.I):
                return "离子转运与稳态"
            return category

    # Phase 1: 关键词匹配（词边界，优先级 first-match-wins）
    kw_cat = _keyword_match(desc)

    # Phase 2: GO 层级辅助
    hier_result = _hierarchy_match(go_id)
    hier_cat = hier_result[0] if hier_result else None
    hier_score = hier_result[1] if hier_result else 0

    # Phase 3: 裁决
    # 层级仅在以下情况覆盖关键词结果：
    #   1. 关键词未命中任何具体类别（None 或 "其他/未归类"）
    #   2. 关键词命中了通用背景兜底类（"发育·其他"），层级指向更具体的非背景类别
    # 层级不会覆盖已由关键词分配的具体背景亚类（如"神经系统发育"）或非背景类别
    if kw_cat and kw_cat not in ("其他/未归类", "发育·其他（去分化背景项）"):
        return kw_cat

    # 关键词未命中具体类别，或命中了通用背景兜底
    if hier_cat and hier_score >= 2:
        # 如果关键词命中了通用背景兜底，层级应指向非背景类别才覆盖
        if kw_cat == "发育·其他（去分化背景项）" and hier_cat in _BACKGROUND_CATEGORIES:
            return kw_cat  # 层级也指向背景，保留关键词结果
        return hier_cat

    # The fallback rules apply only after all specific rules and GO anchors.
    if _ORGANELLE_FALLBACK.search(desc):
        return "细胞器与膜组织"
    if _TRANSPORT_FALLBACK.search(desc):
        return "小分子与跨膜转运"
    if _METABOLISM_FALLBACK.search(desc):
        return "代谢（通用/其他）"
    if _RESPONSE_FALLBACK.search(desc):
        return "刺激与激素响应（其他）"
    return kw_cat or "其他/未归类"


def assign_fine_category(desc, go_id=None):
    """
    将一条 GO 描述归入英文功能类别。

    中文类别可用 assign_fine_category_cn(desc, go_id) 获取。
    """
    return _category_en(assign_fine_category_cn(desc, go_id))


def identify_biological_theme(desc):
    """将一条 GO 描述归入 5 个英文主题之一，不匹配则返回 None。"""
    d = re.sub(r"^\s*obsolete\s+", "", str(desc), flags=re.I).lower()
    for theme, kws in THEME_RULES:
        for kw in kws:
            if kw in d:
                return theme
    return None


def is_signaling_category(category):
    """判断英文类别是否属于指定信号转导类别。"""
    return category in _SIGNALING_CATEGORIES


def assign_merged_category(category, description, background=False):
    """Map a fine category to one merged category without changing Category."""
    if background and category != "Wound Healing & Regeneration":
        return "Development, neural and organogenesis"

    if category in _DNA_GENOME_STABILITY_CATEGORIES:
        return "DNA Replication, Repair & Genome Stability"

    if category in _MIGRATION_ADHESION_CYTOSKELETON_CATEGORIES:
        return "Chemotaxis, Migration, Adhesion & Cytoskeleton"

    if category in _GENERAL_METABOLISM_CATEGORIES:
        return "General Metabolism"

    return category

CLASSIFICATION_COLUMNS = (
    "Category",
    "Category_CN",
    "Background",
    "Signaling",
    "Category_Merged",
    "Theme",
)


def classify_go_bp_term(description, go_id=None):
    """Classify one GO BP term; GO ancestry is supplementary to name rules."""
    go_id = None if pd.isna(go_id) else str(go_id).strip() or None
    description = "" if pd.isna(description) else str(description)
    category_cn = assign_fine_category_cn(description, go_id)
    category = _category_en(category_cn)
    background = category_cn in _BACKGROUND_CATEGORIES
    return {
        "Category": category,
        "Category_CN": category_cn,
        "Background": background,
        "Signaling": is_signaling_category(category),
        "Category_Merged": assign_merged_category(category, description, background),
        "Theme": identify_biological_theme(description),
    }


def classify_go_bp_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of a table with GO BP classification columns added.

    Only ID and Description are required. Other columns and row order are
    preserved; existing classification columns are recalculated.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Expected a pandas DataFrame")

    missing = [col for col in ("ID", "Description") if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    classified = df.copy()
    results = [
        classify_go_bp_term(description, go_id)
        for go_id, description in df[["ID", "Description"]].itertuples(index=False, name=None)
    ]
    for col in CLASSIFICATION_COLUMNS:
        classified[col] = [row[col] for row in results]
    classified[["Background", "Signaling"]] = classified[["Background", "Signaling"]].astype(bool)
    return classified


# Existing callers can migrate to the descriptive names without changing outputs.
categorize_cn = assign_fine_category_cn
categorize = assign_fine_category
theme_of = identify_biological_theme
is_signaling = is_signaling_category
merged_category = assign_merged_category
classify_go_term = classify_go_bp_term
classify_go_terms = classify_go_bp_table
