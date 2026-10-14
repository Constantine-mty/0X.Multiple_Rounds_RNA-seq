"""
GO 本体层级解析模块
====================
解析 go-basic.obo 文件，提取 biological_process 命名空间的术语及其 is_a 父子关系，
计算每个 GO ID 的祖先集合，并缓存到 pickle 文件以加速后续调用。

用法::
    from gobp_pipeline.go_ontology import GOHierarchy
    gh = GOHierarchy()           # 自动下载/加载缓存
    ancestors = gh.get_ancestors("GO:0043039")
    anchor_cats = gh.match_anchors("GO:0043039")  # 返回该术语祖先命中的类别锚点

依赖：无外部依赖（仅标准库 + 可选网络下载）
"""
import os
import re
import pickle
from collections import deque

# go-basic.obo 下载地址
OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
# Keep using the existing project cache after moving this module into the package.
# CACHE_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "go_bp_ancestors.pkl")
CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "go_bp_ancestors.pkl")

# ============================================================================
# 类别锚点 GO ID —— 高层级 GO 术语，若某术语的祖先集合包含这些 ID，
# 则强烈暗示该术语属于对应类别。
# 所有锚点均已验证为当前 GO 本体中存在的、非 obsolete 的术语。
# weight: 锚点置信权重（越高越强）
# ============================================================================
ANCHORS = {
    "翻译与核糖体生物发生": [
        ("GO:0006412", 3),   # translation
        ("GO:0042254", 3),   # ribosome biogenesis
        ("GO:0002181", 2),   # cytoplasmic translation
        ("GO:0070125", 2),   # mitochondrial translational elongation
        ("GO:0006414", 2),   # translational elongation
        ("GO:0006413", 2),   # translational initiation
        ("GO:0006415", 2),   # translational termination
        ("GO:0006417", 2),   # regulation of translation
        ("GO:0043039", 3),   # tRNA aminoacylation（GO层级中不在translation下，但生物学上属于翻译）
    ],
    "RNA加工与剪接": [
        ("GO:0006396", 3),   # RNA processing
        ("GO:0008380", 3),   # RNA splicing
        ("GO:0000375", 2),   # RNA splicing, via transesterification
        ("GO:0006397", 2),   # mRNA processing
        ("GO:0006400", 2),   # tRNA modification
        ("GO:0000154", 2),   # rRNA modification
    ],
    "DNA复制": [
        ("GO:0006260", 3),   # DNA replication
        ("GO:0006270", 2),   # DNA replication initiation
    ],
    "DNA损伤应答与修复": [
        ("GO:0006281", 3),   # DNA repair
        ("GO:0006974", 3),   # cellular response to DNA damage stimulus
        ("GO:0000724", 2),   # double-strand break repair via homologous recombination
        ("GO:0000725", 2),   # recombinational repair
        ("GO:0006289", 2),   # nucleotide-excision repair
        ("GO:0006298", 2),   # mismatch repair
    ],
    "有丝分裂与染色体分离": [
        ("GO:0000278", 3),   # mitotic cell cycle
        ("GO:0140014", 3),   # mitotic nuclear division
        ("GO:0140013", 2),   # meiotic nuclear division
        ("GO:0007059", 2),   # chromosome segregation
        ("GO:0000819", 2),   # sister chromatid segregation
    ],
    "转录与染色质": [
        ("GO:0006351", 3),   # DNA-templated transcription
        ("GO:0006355", 3),   # regulation of DNA-templated transcription
        ("GO:0006325", 3),   # chromatin organization
        ("GO:0006338", 2),   # chromatin remodeling
    ],
    "细胞周期调控与增殖": [
        ("GO:0008283", 3),   # cell population proliferation
        ("GO:0022402", 3),   # cell cycle process
        ("GO:0022403", 2),   # cell cycle phase
        ("GO:0008284", 2),   # positive regulation of cell population proliferation
    ],
    "凋亡、程序性死亡与自噬": [
        ("GO:0006915", 3),   # apoptotic process
        ("GO:0042981", 2),   # regulation of apoptotic process
        ("GO:0008219", 3),   # cell death
        ("GO:0006914", 3),   # autophagy
        ("GO:0016236", 2),   # macroautophagy
        ("GO:0070266", 2),   # necroptotic process
    ],
    "T细胞激活与信号": [
        ("GO:0042110", 3),   # T cell activation
        ("GO:0046634", 2),   # regulation of alpha-beta T cell activation
        ("GO:0046632", 2),   # alpha-beta T cell differentiation
        ("GO:0043367", 2),   # CD4-positive, alpha-beta T cell differentiation
        ("GO:0030217", 2),   # T cell differentiation
    ],
    "B细胞与抗体": [
        ("GO:0042113", 3),   # B cell activation
        ("GO:0030183", 2),   # B cell differentiation
    ],
    "白细胞激活与迁移": [
        ("GO:0045321", 3),   # leukocyte activation
        ("GO:0002283", 2),   # neutrophil activation involved in immune response
        ("GO:0030593", 2),   # neutrophil chemotaxis
        ("GO:0006935", 2),   # chemotaxis
    ],
    "杀伤功能": [
        ("GO:0001909", 3),   # leukocyte mediated cytotoxicity
    ],
    "免疫效应与耐受": [
        ("GO:0006952", 3),   # defense response
        ("GO:0002376", 3),   # immune system process
        ("GO:0006955", 3),   # immune response
        ("GO:0002507", 2),   # tolerance induction
    ],
    "细胞因子与白细胞介素响应": [
        ("GO:0034097", 3),   # response to cytokine
        ("GO:0001817", 2),   # regulation of cytokine production
        ("GO:0001816", 2),   # cytokine production
    ],
    "抗原加工与递呈": [
        ("GO:0019882", 3),   # antigen processing and presentation
        ("GO:0048002", 2),   # antigen processing and presentation of peptide antigen
    ],
    "应激反应": [
        ("GO:0006950", 3),   # response to stress
        ("GO:0009617", 2),   # response to bacterium
    ],
    "干扰素与抗病毒反应": [
        ("GO:0009615", 3),   # response to virus
        ("GO:0051607", 2),   # defense response to virus
        ("GO:0035455", 2),   # response to interferon-alpha
        ("GO:0035456", 2),   # response to interferon-beta
        ("GO:0035458", 2),   # cellular response to interferon-beta
    ],
    # --- 代谢 ---
    "脂质与固醇代谢": [
        ("GO:0006629", 3),   # lipid metabolic process
        ("GO:0008203", 2),   # cholesterol metabolic process
        ("GO:0006633", 2),   # fatty acid biosynthetic process
    ],
    "糖代谢与糖酵解": [
        ("GO:0006006", 3),   # glucose metabolic process
        ("GO:0006096", 3),   # glycolytic process
        ("GO:0005975", 2),   # carbohydrate metabolic process
    ],
    "氨基酸代谢": [
        ("GO:0006520", 3),   # cellular amino acid metabolic process
        ("GO:0009063", 2),   # amino acid catabolic process
    ],
    "核苷酸与嘌呤代谢": [
        ("GO:0009117", 3),   # nucleotide metabolic process
        ("GO:0006163", 2),   # purine nucleotide metabolic process
        ("GO:0006221", 2),   # pyrimidine nucleotide biosynthetic process
    ],
    # --- 背景亚类 ---
    "神经系统发育": [
        ("GO:0007399", 3),   # nervous system development
        ("GO:0022008", 3),   # neurogenesis
        ("GO:0048666", 2),   # neuron development
        ("GO:0007409", 2),   # axonogenesis
        ("GO:0042552", 2),   # myelination
        ("GO:0042063", 2),   # gliogenesis
        ("GO:0007268", 2),   # chemical synaptic transmission
        ("GO:0050808", 2),   # synapse organization
    ],
    "心脏发育": [
        ("GO:0007507", 3),   # heart development
    ],
    "血管与内皮发育": [
        ("GO:0001568", 3),   # blood vessel development
        ("GO:0001525", 3),   # angiogenesis
        ("GO:0001570", 2),   # vasculogenesis
        ("GO:0003158", 2),   # endothelium development
        ("GO:0048514", 2),   # blood vessel morphogenesis
    ],
    "骨骼与软骨发育": [
        ("GO:0001501", 3),   # skeletal system development
        ("GO:0001503", 2),   # ossification
        ("GO:0051216", 2),   # cartilage development
        ("GO:0002062", 2),   # chondrocyte differentiation
    ],
    "肌肉发育": [
        ("GO:0007517", 3),   # muscle organ development
        ("GO:0042692", 2),   # muscle cell differentiation
        ("GO:0014706", 2),   # striated muscle tissue development
    ],
    "上皮与皮肤发育": [
        ("GO:0008544", 3),   # epidermis development
        ("GO:0048730", 2),   # epidermis morphogenesis
        ("GO:0060429", 2),   # epithelium development
        ("GO:0002064", 2),   # epithelial cell development
    ],
    "器官发育": [
        ("GO:0048513", 3),   # animal organ development
        ("GO:0001822", 2),   # kidney development
        ("GO:0001654", 2),   # eye development
        ("GO:0048568", 2),   # embryonic organ development
    ],
    "胚胎与体轴发育": [
        ("GO:0009790", 3),   # embryo development
        ("GO:0001704", 2),   # formation of primary germ layer
        ("GO:0007369", 2),   # gastrulation
        ("GO:0009950", 2),   # dorsal/ventral axis specification
        ("GO:0009952", 2),   # anterior/posterior pattern specification
    ],
    "形态发生": [
        ("GO:0009653", 3),   # anatomical structure morphogenesis
        ("GO:0032989", 2),   # cellular anatomical entity morphogenesis
        ("GO:0035239", 2),   # tube morphogenesis
    ],
    "细胞分化与命运决定": [
        ("GO:0030154", 3),   # cell differentiation
        ("GO:0001708", 2),   # cell fate specification
        ("GO:0001709", 2),   # cell fate determination
        ("GO:0045165", 2),   # cell fate commitment
    ],
    "伤口愈合与再生": [
        ("GO:0042060", 3),   # wound healing
        ("GO:0031099", 2),   # regeneration
        ("GO:0048771", 2),   # tissue remodeling
    ],
    "系统生理与行为": [
        ("GO:0007585", 2),   # respiratory gaseous exchange by respiratory system
        ("GO:0003093", 2),   # regulation of glomerular filtration
        ("GO:0008015", 2),   # blood circulation
        ("GO:0007589", 2),   # body fluid secretion
        ("GO:0007610", 2),   # behavior
        ("GO:0050817", 2),   # coagulation
        ("GO:0007599", 2),   # hemostasis
        ("GO:0044057", 2),   # regulation of system process
    ],
}


class GOHierarchy:
    """GO 本体层级解析与祖先查询。"""

    def __init__(self, obo_path=None, cache_path=None, force_download=False):
        self.cache_path = cache_path or CACHE_PATH
        self.obo_path = obo_path
        self.terms = {}
        self.ancestors = {}
        self.id_to_name = {}
        self.name_to_id = {}
        self._loaded = False
        self._load(force_download)

    def _load(self, force_download=False):
        # 尝试从缓存加载
        if not force_download and os.path.exists(self.cache_path):
            try:
                with open(self.cache_path, "rb") as f:
                    data = pickle.load(f)
                self.terms = data["terms"]
                self.ancestors = data["ancestors"]
                self.id_to_name = data["id_to_name"]
                self.name_to_id = data["name_to_id"]
                self._loaded = True
                return
            except Exception:
                pass

        # 确定 OBO 文件路径
        obo_path = self.obo_path
        if obo_path is None or not os.path.exists(obo_path):
            obo_path = os.path.join(os.path.dirname(self.cache_path), "go-basic.obo")

        if not os.path.exists(obo_path) or force_download:
            self._download_obo(obo_path)

        self._parse_obo(obo_path)
        self._compute_ancestors()
        self.id_to_name = {gid: t["name"] for gid, t in self.terms.items()}
        self.name_to_id = {t["name"].lower(): gid for gid, t in self.terms.items()}

        try:
            with open(self.cache_path, "wb") as f:
                pickle.dump({
                    "terms": self.terms,
                    "ancestors": self.ancestors,
                    "id_to_name": self.id_to_name,
                    "name_to_id": self.name_to_id,
                }, f, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            pass

        self._loaded = True

    def _download_obo(self, dest_path):
        import urllib.request
        urllib.request.urlretrieve(OBO_URL, dest_path)

    def _parse_obo(self, obo_path):
        terms = {}
        current = None
        with open(obo_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line == "[Term]":
                    if current:
                        terms[current["id"]] = current
                    current = {"id": "", "name": "", "namespace": "", "is_a": []}
                elif current is not None:
                    if line.startswith("id:"):
                        current["id"] = line[4:].strip()
                    elif line.startswith("name:"):
                        current["name"] = line[6:].strip()
                    elif line.startswith("namespace:"):
                        current["namespace"] = line[11:].strip()
                    elif line.startswith("is_a:"):
                        m = re.match(r"is_a:\s+(GO:\d+)", line)
                        if m:
                            current["is_a"].append(m.group(1))
            if current:
                terms[current["id"]] = current
        self.terms = terms

    def _compute_ancestors(self):
        ancestors = {}
        def _get_anc(gid):
            if gid in ancestors:
                return ancestors[gid]
            result = set()
            queue = deque(self.terms.get(gid, {}).get("is_a", []))
            while queue:
                aid = queue.popleft()
                if aid in result:
                    continue
                result.add(aid)
                queue.extend(self.terms.get(aid, {}).get("is_a", []))
            ancestors[gid] = result
            return result
        for gid in self.terms:
            _get_anc(gid)
        self.ancestors = ancestors

    def get_ancestors(self, go_id):
        return self.ancestors.get(go_id, set())

    def get_ancestor_names(self, go_id):
        ancs = self.get_ancestors(go_id)
        return [self.id_to_name[a] for a in ancs if a in self.id_to_name]

    def match_anchors(self, go_id):
        """检查某 GO ID 的祖先集合是否命中类别锚点。返回 {category: total_weight}。"""
        if not go_id or go_id not in self.ancestors:
            return {}
        ancs = self.ancestors[go_id]
        all_ids = ancs | {go_id}
        scores = {}
        for cat, anchors in ANCHORS.items():
            total = 0
            for aid, weight in anchors:
                if aid in all_ids:
                    total += weight
            if total > 0:
                scores[cat] = total
        return scores

    def is_descendant_of(self, go_id, ancestor_go_id):
        if not go_id:
            return False
        return ancestor_go_id in self.get_ancestors(go_id)

    def lookup_id(self, description):
        return self.name_to_id.get(description.lower().strip())


if __name__ == "__main__":
    gh = GOHierarchy()
    print(f"Loaded {len(gh.terms)} terms, {len(gh.ancestors)} ancestor sets")

    # 验证所有锚点存在且非 obsolete
    print("\n=== Anchor verification ===")
    for cat, anchors in ANCHORS.items():
        for aid, weight in anchors:
            name = gh.id_to_name.get(aid, "NOT FOUND")
            obsolete = "obsolete" in name.lower() if name != "NOT FOUND" else True
            flag = " *** OBSOLETE ***" if obsolete else ""
            if obsolete or name == "NOT FOUND":
                print(f"  {cat:25s} {aid} w={weight}  {name}{flag}")

    # 测试
    test_cases = [
        ("GO:0043039", "tRNA aminoacylation"),
        ("GO:0008033", "tRNA processing"),
        ("GO:0000460", "maturation of 5.8S rRNA"),
        ("GO:0014032", "neural crest cell development"),
        ("GO:0045576", "mast cell activation"),
        ("GO:0045444", "fat cell differentiation"),
        ("GO:0045773", "positive regulation of axon extension"),
    ]
    print("\n=== Test cases ===")
    for gid, name in test_cases:
        scores = gh.match_anchors(gid)
        top = sorted(scores.items(), key=lambda x: -x[1])[:3]
        print(f"\n{gid} {name}:")
        for cat, score in top:
            print(f"  {cat}: {score}")
