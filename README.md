---

````markdown
# 🧬 GenBank Gene Comparison and FASTA Extractor

This Python script compares **genes across multiple GenBank files**, identifies shared genes among all species, and generates:

- A detailed **HTML report** summarizing gene presence across files
- Individual **FASTA files** for each shared gene, ready for alignment

Perfect for comparative genomics, phylogenetics, or chloroplast genome analysis.

---

## 📌 Features

- Accepts multiple GenBank files or folders
- Extracts gene information: name, locus tag, product, and location
- Identifies genes shared across all species
- Outputs FASTA files for each shared gene
- Generates an interactive HTML report with:
  - Presence/absence matrix
  - Summary statistics
  - Gene metadata per file

---

## 🛠️ Requirements

- Python 3.6+
- Biopython

Install dependencies:

```bash
pip install biopython
````

---

## 🚀 Usage

```bash
python compare_genbank_genes.py PATH1 [PATH2 ...] --reference REFERENCE_FILE [--output OUTPUT.html]
```

### Arguments:

* `PATH1 [PATH2 ...]`: One or more GenBank files or directories containing `.gb`, `.gbk`, `.genbank`, or `.gbf` files
* `--reference`: Reference GenBank file (filename only, not path) to compare others against
* `--output`: (Optional) Output HTML file name (default: `genbank_comparison_report.html`)

---

### 🔍 Example:

```bash
python compare_genbank_genes.py ./genomes/ --reference NC_030785.gbk
```

This will:

* Search all valid GenBank files in `./genomes/`
* Compare genes to `NC_030785.gbk`
* Identify shared genes
* Save FASTA files for each shared gene in `./shared_genes_fasta/`
* Generate `genbank_comparison_report.html`

---

## 📂 Output

* `genbank_comparison_report.html`: Interactive HTML report
* `shared_genes_fasta/`: Directory with one FASTA file per shared gene

Example FASTA header format:

```text
>NC_030785_1_petG cytochrome b6/f complex subunit V _ Artemisia argyi chloroplast complete genome
```

---

## 🧪 Applications

* Comparative genomics
* Organellar genome analysis
* Phylogenetic studies
* Gene conservation and annotation review

---

## 📄 License

MIT License

---

## 🙋‍♂️ Author

Developed by \[Aymen Trso]. Contributions and feedback are welcome!

```

---
