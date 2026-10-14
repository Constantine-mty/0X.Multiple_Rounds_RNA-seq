from pathlib import Path

import pandas as pd
from gobp_category import go_bp_classification as bp
from gobp_category.gsea_classification import classify_gsea_results

if not bp._HIERARCHY_AVAILABLE:
    raise RuntimeError("GO 层级未加载，请检查 go_bp_ancestors.pkl 的路径")


from pathlib import Path
import pandas as pd
from gobp_category.gsea_classification import classify_gsea_results

source = Path("../table/enrichment_results/gsea_sig_group.merged.csv")
output = source.with_name("gsea_sig_group.merged.classified.csv")

df = pd.read_csv(source, encoding="utf-8-sig")
if "ID" not in df.columns:
    raise ValueError("缺少 ID 列；没有 GO ID 就无法利用层级关系")
classified = classify_gsea_results(df)
classified.to_csv(output, index=False, encoding="utf-8-sig")
print(f"已保存 {len(classified)} 行：{output}")






