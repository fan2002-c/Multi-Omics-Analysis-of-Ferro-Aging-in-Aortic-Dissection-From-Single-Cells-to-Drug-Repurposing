# GitHub 仓库上传指南

## 状态：本地仓库完全就绪，等待推送

所有分析代码已整理完毕，Git 仓库已在本地创建并提交。
由于工作环境网络限制，无法直接推送到 GitHub，请按以下步骤操作。

---

## 方式一：手把手操作（推荐）

### 1. 创建 GitHub 仓库
打开 https://github.com/new
- Repository name: `ferro-aging-ad`
- Description: `Multi-omics characterization of ferro-aging in aortic dissection: scRNA-seq, bulk validation, drug repurposing, and cross-species age analysis`
- 选择 **Public**（学术期刊要求公开）
- **不要**勾选 "Add a README file"（我们已有 README）
- 点击 "Create repository"

### 2. 推送代码
打开 PowerShell 或 Git Bash，执行：

```bash
# 进入仓库目录
cd "C:\Users\lidaf\WorkBuddy\2026-06-11-19-29-35\github_repo"

# 添加远程仓库（替换 YOUR_USERNAME）
git remote add origin https://github.com/YOUR_USERNAME/ferro-aging-ad.git

# 推送
git push -u origin main
```

### 3. 更新手稿中的链接
推送成功后，将 `submission_package/manuscript.md` 中的
`https://github.com/lidaf/ferro-aging-ad` 替换为实际仓库 URL。

---

## 方式二：上传 ZIP（备选）

如果无法使用 Git 命令行：
1. 打开 https://github.com/new 创建仓库（同上）
2. 在仓库页面点击 "uploading an existing file"
3. 将 `submission_package/ferro_aging_code_repo.zip` 拖入上传区域
4. 提交

---

## 仓库内容清单

```
ferro-aging-ad/
├── README.md           → 项目说明（英文，投稿级）
├── LICENSE             → MIT License
├── requirements.txt    → Python 依赖
├── .gitignore          → Git 忽略规则
├── data/
│   └── ferro_aging_geneset.csv    → 129基因 ferro-aging 基因集
├── scripts/
│   ├── build_ferro_aging_geneset.py     → 01 基因集构建
│   ├── ferro_aging_main_analysis.py     → 02 scRNA-seq 主分析
│   ├── ferro_aging_deep_analysis.py     → 03 亚群深度分析
│   ├── ferro_aging_bulk_validation.py   → 04 Bulk RNA-seq 验证
│   ├── dgidb_query.py                   → 05 DGIdb 药物查询
│   ├── ferro_aging_drug_prediction.py   → 06 药物优先级排序
│   ├── parse_primate_aging.py           → 07a 灵长类数据解析
│   ├── ferro_aging_age_validation.py    → 07b 年龄交叉验证
│   └── manuscript_figures.py            → 08 出版级图表生成
└── results/            → 空目录，分析结果输出位置
```

---

## 文件位置

- **Git 仓库（本地）**：`github_repo/`
- **ZIP 包**：`submission_package/ferro_aging_code_repo.zip`
- **手稿（已更新 Code Availability）**：`submission_package/manuscript.md`
