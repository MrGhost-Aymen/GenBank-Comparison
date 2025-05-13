import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os
from datetime import datetime

def find_genbank_files(input_paths):
    """Find all GenBank files from input paths (both files and directories)"""
    genbank_files = []
    valid_extensions = ('.gb', '.gbk', '.genbank', '.gbf')
    
    for path in input_paths:
        if os.path.isfile(path):
            if path.lower().endswith(valid_extensions):
                genbank_files.append(path)
            else:
                print(f"Warning: {path} is not a GenBank file (skipping)")
        elif os.path.isdir(path):
            for root, _, files in os.walk(path):
                for file in files:
                    if file.lower().endswith(valid_extensions):
                        genbank_files.append(os.path.join(root, file))
        else:
            print(f"Warning: {path} is neither a file nor directory (skipping)")
    
    if not genbank_files:
        raise ValueError(f"No valid GenBank files found in input paths: {', '.join(input_paths)}")
    
    return genbank_files

def parse_genbank_files(file_paths):
    gene_data = {}
    all_sequences = defaultdict(dict)
    
    for file_path in file_paths:
        try:
            records = list(SeqIO.parse(file_path, "genbank"))
            file_name = os.path.basename(file_path)
            file_genes = set()
            file_gene_details = []
            
            for record in records:
                organism = record.annotations.get('organism', 'Unknown')
                molecule_type = record.annotations.get('molecule_type', 'Unknown')
                
                for feature in record.features:
                    if feature.type == "gene":
                        gene_name = feature.qualifiers.get('gene', ['unnamed'])[0]
                        locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                        product = feature.qualifiers.get('product', [''])[0]
                        file_genes.add(gene_name)
                        
                        # Store gene details
                        file_gene_details.append({
                            'name': gene_name,
                            'locus_tag': locus_tag,
                            'product': product,
                            'location': str(feature.location)
                        })
                        
                        # Extract and store sequence
                        try:
                            sequence = str(feature.extract(record.seq))
                            if sequence:
                                all_sequences[gene_name][file_name] = {
                                    'sequence': sequence,
                                    'locus_tag': locus_tag,
                                    'product': product,
                                    'organism': organism,
                                    'molecule_type': molecule_type,
                                    'accession': record.id,
                                    'definition': record.description
                                }
                        except Exception as e:
                            print(f"Warning: Could not extract sequence for gene {gene_name} in {file_name}: {str(e)}")
            
            gene_data[file_name] = {
                'genes': file_genes,
                'gene_details': file_gene_details
            }
        except Exception as e:
            print(f"Warning: Could not parse {file_path} (skipping). Error: {str(e)}")
    
    return gene_data, all_sequences

def find_shared_genes(gene_data):
    """Find genes present in all species"""
    if not gene_data:
        return set()
    
    # Start with genes from first file
    shared_genes = set(gene_data[next(iter(gene_data))]['genes'])
    
    # Intersect with genes from other files
    for file_name, data in gene_data.items():
        shared_genes.intersection_update(data['genes'])
    
    return shared_genes

def save_shared_genes_fasta(shared_genes, all_sequences, output_dir='shared_genes_fasta'):
    """
    Save shared genes in FASTA format with specific header format:
    ">174369758543892_petG cytochrome b6/f complex subunit V | Artemisia ifranensis chloroplast, complete genome"
    One FASTA file per gene named like "petG.fasta"
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for gene in shared_genes:
        gene_sequences = all_sequences.get(gene, {})
        if not gene_sequences:
            print(f"Warning: No sequences found for shared gene {gene}")
            continue
        
        # Create filename using only the gene name
        fasta_filename = os.path.join(output_dir, f"{gene}.fasta")
        records = []
        
        for file_name, data in gene_sequences.items():
            # Create header in the specified format
            header = f"{data['accession']}_{data['locus_tag']}_{gene} {data['product']} | {data['definition']}"
            
            # Clean up header to remove problematic characters
            header = header.replace(",", "").replace(":", "").replace(";", "").replace("|", "_").strip()
            
            # Create SeqRecord
            seq_record = SeqRecord(
                Seq(data['sequence']),
                id=header,
                description=""
            )
            records.append(seq_record)
        
        # Write all sequences for this gene to one FASTA file
        with open(fasta_filename, 'w') as handle:
            SeqIO.write(records, handle, 'fasta')
        
        print(f"Saved {len(records)} sequences for gene {gene} to {fasta_filename}")
    
    print(f"\nAll shared genes saved in FASTA format to '{output_dir}' directory")
    print(f"Each file contains sequences from all species for one gene, ready for alignment")

def compare_genes(gene_data, reference_file):
    """Compare genes across files with respect to a reference file"""
    if reference_file not in gene_data:
        raise ValueError(f"Reference file {reference_file} not found in input files")
    
    reference_genes = gene_data[reference_file]['genes']
    all_genes = set(reference_genes)
    
    # Collect all genes from all files
    for file_name, data in gene_data.items():
        all_genes.update(data['genes'])
    
    # Create comparison matrix
    comparison = []
    for gene in sorted(all_genes):
        gene_info = {'gene': gene}
        for file_name in gene_data.keys():
            gene_info[file_name] = gene in gene_data[file_name]['genes']
        comparison.append(gene_info)
    
    return comparison

def generate_html_report(gene_data, comparison, reference_file, output_file, shared_genes):
    """Generate an HTML report from the comparison data"""
    file_names = list(gene_data.keys())
    
    # Get current date for report
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Use HTML entities instead of Unicode characters
    check_mark = "✔"  # ✓
    cross_mark = "✘"  # ✗
    
    # Start HTML with explicit UTF-8 encoding
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>GenBank File Comparison Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                margin: 20px;
                color: #333;
            }}
            h1, h2 {{
                color: #2c3e50;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
                position: sticky;
                top: 0;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .present {{
                background-color: #d4edda;
                text-align: center;
            }}
            .absent {{
                background-color: #f8d7da;
                text-align: center;
            }}
            .reference {{
                font-weight: bold;
                background-color: #fff3cd;
            }}
            .shared {{
                font-weight: bold;
                background-color: #cce5ff;
            }}
            .summary {{
                margin-bottom: 30px;
                padding: 15px;
                background-color: #e2e3e5;
                border-radius: 5px;
            }}
            .gene-details {{
                margin-top: 30px;
            }}
            .section {{
                margin-top: 40px;
            }}
            .fasta-info {{
                background-color: #e8f4f8;
                padding: 15px;
                border-radius: 5px;
                margin: 20px 0;
            }}
        </style>
    </head>
    <body>
        <h1>GenBank File Comparison Report</h1>
        <p>Generated on: {current_date}</p>
        <p>Reference file: <strong>{reference_file}</strong></p>
        
        <div class="summary">
            <h2>Summary</h2>
            <p>Number of files compared: {len(file_names)}</p>
            <p>Total unique genes found: {len(comparison)}</p>
            <p>Genes in reference: {len(gene_data[reference_file]['genes'])}</p>
            <p>Genes shared by all species: {len(shared_genes)}</p>
        </div>
        
        <div class="fasta-info">
            <h2>FASTA Files for Alignment</h2>
            <p>Sequences for shared genes have been saved in FASTA format in the 'shared_genes_fasta' directory.</p>
            <p>Each gene has its own FASTA file (e.g., 'petG.fasta') containing sequences from all species with headers formatted as:</p>
            <p><code>>NC_030785_1_petG cytochrome b6/f complex subunit V | Artemisia argyi chloroplast, complete genome</code></p>
            <p>These files are ready for alignment with tools like MUSCLE, MAFFT, or Clustal.</p>
        </div>
    """
    
    # Add shared genes section
    html += """
        <div class="section">
            <h2>Genes Shared by All Species</h2>
            <ul>
    """
    for gene in sorted(shared_genes):
        html += f'<li>{gene}</li>'
    html += """
            </ul>
        </div>
    """
    
    # Add comparison table
    html += """
        <div class="section">
            <h2>Gene Presence Comparison</h2>
            <table>
                <thead>
                    <tr>
                        <th>Gene</th>
    """
    for file_name in file_names:
        if file_name == reference_file:
            html += f'<th class="reference">{file_name} (reference)</th>'
        else:
            html += f'<th>{file_name}</th>'
    html += """
                    </tr>
                </thead>
                <tbody>
    """
    
    for gene_info in comparison:
        is_shared = gene_info["gene"] in shared_genes
        row_class = "shared" if is_shared else ""
        html += f'<tr class="{row_class}"><td>{gene_info["gene"]}</td>'
        for file_name in file_names:
            status = "present" if gene_info[file_name] else "absent"
            if file_name == reference_file and not gene_info[file_name]:
                status = "absent reference"
            symbol = check_mark if gene_info[file_name] else cross_mark
            html += f'<td class="{status}">{symbol}</td>'
        html += '</tr>'
    
    html += """
                </tbody>
            </table>
        </div>
    """
    
    # Add detailed gene information for each file
    html += """
        <div class="section gene-details">
            <h2>Detailed Gene Information by File</h2>
    """
    
    for file_name in file_names:
        html += f"""
            <h3>{file_name}</h3>
            <table>
                <thead>
                    <tr>
                        <th>Gene Name</th>
                        <th>Locus Tag</th>
                        <th>Product</th>
                        <th>Location</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for gene in sorted(gene_data[file_name]['gene_details'], key=lambda x: x['name']):
            html += f"""
                <tr>
                    <td>{gene['name']}</td>
                    <td>{gene['locus_tag']}</td>
                    <td>{gene['product']}</td>
                    <td>{gene['location']}</td>
                </tr>
            """
        
        html += """
                </tbody>
            </table>
        """
    
    html += """
        </div>
    </body>
    </html>
    """
    
    # Write to file with explicit UTF-8 encoding
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"HTML report generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Compare GenBank files and generate an HTML report.')
    parser.add_argument('paths', metavar='PATH', type=str, nargs='+',
                        help='GenBank files or directories containing GenBank files to compare')
    parser.add_argument('--reference', type=str, required=True,
                        help='Reference file name (must match one of the input file names)')
    parser.add_argument('--output', type=str, default='genbank_comparison_report.html',
                        help='Output HTML file name (default: genbank_comparison_report.html)')
    
    args = parser.parse_args()
    
    # Find all GenBank files in input paths
    genbank_files = find_genbank_files(args.paths)
    print(f"Found {len(genbank_files)} GenBank files to compare")
    
    # Parse GenBank files
    gene_data, all_sequences = parse_genbank_files(genbank_files)
    
    # Verify reference file exists in the parsed data
    if args.reference not in gene_data:
        # Try to match by basename if full path was given
        reference_basename = os.path.basename(args.reference)
        if reference_basename in gene_data:
            args.reference = reference_basename
        else:
            raise ValueError(f"Reference file {args.reference} not found in input files. "
                           f"Available files: {', '.join(gene_data.keys())}")
    
    # Find shared genes
    shared_genes = find_shared_genes(gene_data)
    print(f"Found {len(shared_genes)} genes shared by all species")
    
    # Save shared genes in FASTA format
    if shared_genes:
        save_shared_genes_fasta(shared_genes, all_sequences)
    
    # Compare genes
    comparison = compare_genes(gene_data, args.reference)
    
    # Generate HTML report
    generate_html_report(gene_data, comparison, args.reference, args.output, shared_genes)

if __name__ == "__main__":
    main()