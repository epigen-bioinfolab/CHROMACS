import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import tempfile
import subprocess
import shutil
import sys
import argparse
from pathlib import Path

class ATACSeqAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("ChromAcS+ Additional Analysis Toolkit")
        self.root.geometry("800x600")
        self.root.minsize(800, 600)
        
        self._make_title_bar()

        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        self.notebook = ttk.Notebook(self.main_frame)
        self.tab_overlap = ttk.Frame(self.notebook)
        self.tab_expression = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_overlap, text="Peak Overlap Analysis")
        self.notebook.add(self.tab_expression, text="Compare with Expression")
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.setup_overlap_tab()
        self.setup_expression_tab()

        self.status = tk.Label(self.main_frame, text="Ready", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

    def _make_title_bar(self):
        BAR_HEIGHT = 40
        self.title_bar = tk.Frame(self.root, bg="maroon", height=BAR_HEIGHT)
        self.title_bar.pack(fill=tk.X)
        self.title_bar.pack_propagate(False)

        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        icon_path = os.path.join(script_dir, "ChromAcS.png")
        if os.path.exists(icon_path):
            try:
                img = tk.PhotoImage(file=icon_path)
                # Scale down if too tall
                max_h = BAR_HEIGHT - 8
                if img.height() > max_h:
                    factor = img.height() // max_h + 1
                    img = img.subsample(factor, factor)
                self.logo_img = img
                logo_label = tk.Label(self.title_bar, image=self.logo_img, bg="maroon")
                logo_label.image = img  # Keep reference
                logo_label.pack(side=tk.LEFT, padx=5, pady=4)
            except Exception as e:
                print(f"Could not load logo: {str(e)}")

        title_font = ("Helvetica", 12, "bold")
        tk.Label(self.title_bar, text=self.root.title(),
                 bg="maroon", fg="white", font=title_font).pack(side=tk.LEFT, padx=5)

    # ------------------------------------------
    # Peak Overlap Analysis Tab
    # ------------------------------------------
    def setup_overlap_tab(self):
        self.tab_overlap.configure(style='Custom.TFrame')

        tk.Label(self.tab_overlap, text="ATAC-seq Peaks File:").grid(row=0, column=0, padx=5, pady=5, sticky='e')
        self.overlap_input = tk.Entry(self.tab_overlap, width=50)
        self.overlap_input.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(self.tab_overlap, text="Browse", command=lambda: self.browse_file(self.overlap_input)).grid(row=0, column=2, padx=5)

        tk.Label(self.tab_overlap, text="Window Size (bp):").grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.window_size = tk.Entry(self.tab_overlap, width=10)
        self.window_size.grid(row=1, column=1, padx=5, pady=5, sticky='w')
        self.window_size.insert(0, "200")

        self.bed_entries = []
        for i in range(5):
            tk.Label(self.tab_overlap, text=f"BED File {i+1}:").grid(row=i+2, column=0, padx=5, pady=5, sticky='e')
            entry = tk.Entry(self.tab_overlap, width=50)
            entry.grid(row=i+2, column=1, padx=5, pady=5)
            tk.Button(self.tab_overlap, text="Browse", command=lambda e=entry: self.browse_file(e)).grid(row=i+2, column=2, padx=5)
            self.bed_entries.append(entry)

        tk.Label(self.tab_overlap, text="Output File:").grid(row=7, column=0, padx=5, pady=5, sticky='e')
        self.overlap_output = tk.Entry(self.tab_overlap, width=50)
        self.overlap_output.grid(row=7, column=1, padx=5, pady=5)
        tk.Button(self.tab_overlap, text="Browse", command=lambda: self.browse_save_file(self.overlap_output)).grid(row=7, column=2, padx=5)

        run_button = tk.Button(self.tab_overlap, text="Run Overlap Analysis", 
                              command=self.run_overlap,
                              bg="#4a7abc", fg="white", 
                              font=("Helvetica", 10, "bold"),
                              padx=10, pady=5)
        run_button.grid(row=8, column=1, pady=15)

    def run_overlap(self):
        input_file = self.overlap_input.get()
        output_file = self.overlap_output.get()
        window = self.window_size.get()
        bed_files = [e.get() for e in self.bed_entries if e.get().strip()]
        if not input_file or not output_file:
            messagebox.showerror("Error", "Please specify input and output files")
            return
        if not window.isdigit():
            messagebox.showerror("Error", "Window size must be a positive integer")
            return
        window = int(window)
        self.update_status("Starting overlap analysis...")
        try:
            temp_dir = tempfile.mkdtemp()
            cleaned_input = os.path.join(temp_dir, "cleaned_input.bed")
            self.clean_input_file(input_file, cleaned_input)
            extended_input = os.path.join(temp_dir, "extended_input.bed")
            self.extend_ranges(cleaned_input, extended_input, window)

            processed_beds = []
            chrom_style = self.get_chrom_style(extended_input)
            for i, bed in enumerate(bed_files):
                outb = os.path.join(temp_dir, f"bed_{i}.bed")
                self.process_bed_file(bed, outb, chrom_style)
                processed_beds.append(outb)

            self.run_bedtools_chain(extended_input, processed_beds, output_file)
            self.update_status(f"Analysis complete! Output saved to {output_file}")
            messagebox.showinfo("Success", "Overlap analysis completed successfully")
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
            self.update_status(f"Error: {str(e)}")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def clean_input_file(self, input_path, output_path):
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                line = line.replace(' ', '_')
                parts = line.strip().split('\t')
                cleaned = [part if part else '----' for part in parts]
                outfile.write('\t'.join(cleaned) + '\n')

    def extend_ranges(self, input_path, output_path, window):
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    chrom = parts[0]
                    start = max(0, int(parts[1]) - window)
                    end = int(parts[2]) + window
                    new_line = [chrom, str(start), str(end)] + parts[3:]
                    outfile.write('\t'.join(new_line) + '\n')
                except ValueError:
                    continue

    def get_chrom_style(self, bed_file):
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                return 'with_chr' if parts and parts[0].startswith('chr') else 'without_chr'

    def process_bed_file(self, input_path, output_path, chrom_style):
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if not parts: continue
                chrom = parts[0]
                if chrom_style == 'with_chr' and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                elif chrom_style == 'without_chr' and chrom.startswith('chr'):
                    chrom = chrom[3:]
                parts[0] = chrom
                outfile.write('\t'.join(parts) + '\n')

    def run_bedtools_chain(self, input_file, bed_files, output_file):
        cmd = f"bedtools intersect -a {input_file} -b {bed_files[0]} -wao"
        for b in bed_files[1:]:
            cmd += f" | bedtools intersect -a stdin -b {b} -wao"
        cmd += f" > {output_file}"
        subprocess.run(cmd, shell=True, check=True)

    # ------------------------------------------
    # Expression Comparison Tab
    # ------------------------------------------
    def setup_expression_tab(self):
        self.tab_expression.configure(style='Custom.TFrame')
        
        tk.Label(self.tab_expression, text="Annotated ATAC-seq File:").grid(row=0, column=0, padx=5, pady=5, sticky='e')
        self.atac_input = tk.Entry(self.tab_expression, width=50)
        self.atac_input.grid(row=0, column=1, padx=5, pady=5)
        tk.Button(self.tab_expression, text="Browse", command=lambda: self.browse_file(self.atac_input)).grid(row=0, column=2, padx=5)

        tk.Label(self.tab_expression, text="Expression Data File:").grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.expr_input = tk.Entry(self.tab_expression, width=50)
        self.expr_input.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(self.tab_expression, text="Browse", command=lambda: self.browse_file(self.expr_input)).grid(row=1, column=2, padx=5)

        tk.Label(self.tab_expression, text="ATAC-seq ID Column:").grid(row=2, column=0, padx=5, pady=5, sticky='e')
        self.atac_col = tk.Entry(self.tab_expression, width=5)
        self.atac_col.grid(row=2, column=1, padx=5, pady=5, sticky='w')
        self.atac_col.insert(0, "17")

        tk.Label(self.tab_expression, text="Expression ID Column:").grid(row=3, column=0, padx=5, pady=5, sticky='e')
        self.expr_col = tk.Entry(self.tab_expression, width=5)
        self.expr_col.grid(row=3, column=1, padx=5, pady=5, sticky='w')
        self.expr_col.insert(0, "1")

        tk.Label(self.tab_expression, text="Output File:").grid(row=4, column=0, padx=5, pady=5, sticky='e')
        self.expr_output = tk.Entry(self.tab_expression, width=50)
        self.expr_output.grid(row=4, column=1, padx=5, pady=5)
        tk.Button(self.tab_expression, text="Browse", command=lambda: self.browse_save_file(self.expr_output)).grid(row=4, column=2, padx=5)

        run_button = tk.Button(self.tab_expression, text="Run Expression Merge", 
                              command=self.run_expression,
                              bg="#4a7abc", fg="white", 
                              font=("Helvetica", 10, "bold"),
                              padx=10, pady=5)
        run_button.grid(row=5, column=1, pady=15)

    def run_expression(self):
        atac_file = self.atac_input.get()
        expr_file = self.expr_input.get()
        output_file = self.expr_output.get()
        atac_col = self.atac_col.get()
        expr_col = self.expr_col.get()

        if not atac_file or not expr_file or not output_file:
            messagebox.showerror("Error", "All file fields are required")
            return
        if not atac_col.isdigit() or not expr_col.isdigit():
            messagebox.showerror("Error", "Column indices must be integers")
            return

        atac_idx = int(atac_col) - 1
        expr_idx = int(expr_col) - 1
        self.update_status("Merging expression data...")

        try:
            self.append_expression_data(atac_file, expr_file, output_file, atac_idx, expr_idx)
            self.update_status(f"Merge complete! Output saved to {output_file}")
            messagebox.showinfo("Success", "Expression data merged successfully")
        except Exception as e:
            messagebox.showerror("Error", f"Merge failed: {str(e)}")
            self.update_status(f"Error: {str(e)}")

    def append_expression_data(self, input_file, append_file, output_file, input_col, append_col, na_val="NA"):
        expr_data = {}
        expr_columns = []

        with open(append_file, 'r') as fin:
            header = fin.readline().strip().split('\t')
            expr_columns = header[append_col + 1:]  # Get all columns after the gene symbol

            for line in fin:
                parts = line.strip().split('\t')
                if len(parts) <= append_col:
                    continue
                gene_symbol = parts[append_col]
                expr_values = parts[append_col + 1:]
                expr_data[gene_symbol] = expr_values

        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            header = fin.readline().strip().split('\t')

            new_header = header[:input_col + 1] + expr_columns + header[input_col + 1:]
            fout.write('\t'.join(new_header) + '\n')

            for line in fin:
                parts = line.strip().split('\t')
                if len(parts) <= input_col:
                    fout.write('\t'.join(parts) + '\n')
                    continue

                gene_symbol = parts[input_col]
                if gene_symbol in expr_data:
                    new_parts = parts[:input_col + 1] + expr_data[gene_symbol] + parts[input_col + 1:]
                else:
                    new_parts = parts[:input_col + 1] + [na_val] * len(expr_columns) + parts[input_col + 1:]

                fout.write('\t'.join(new_parts) + '\n')

    # ------------------------------------------
    # Common GUI Helpers
    # ------------------------------------------
    def browse_file(self, entry):
        path = filedialog.askopenfilename()
        
        if path:
            entry.delete(0, tk.END)
            entry.insert(0, path)

    def browse_save_file(self, entry):
        path = filedialog.asksaveasfilename()
        if path:
            entry.delete(0, tk.END)
            entry.insert(0, path)

    def update_status(self, msg):
        self.status.config(text=msg)
        self.status.update_idletasks()


def cli_overlap(args):
    input_file = args.input
    output_file = args.output
    window = args.window
    bed_files = args.beds

    if not input_file or not output_file or not bed_files:
        print("Need --input, --output, and at least one --bed")
        return 1
    if window < 0:
        print("Window must be non-negative")
        return 1

    analyzer = ATACSeqAnalyzer.__new__(ATACSeqAnalyzer)
    analyzer.status = type("S", (), {"config": lambda *a, **k: None, "update_idletasks": lambda: None})()

    temp_dir = tempfile.mkdtemp()
    try:
        cleaned_input = os.path.join(temp_dir, "cleaned_input.bed")
        analyzer.clean_input_file(input_file, cleaned_input)
        extended_input = os.path.join(temp_dir, "extended_input.bed")
        analyzer.extend_ranges(cleaned_input, extended_input, window)

        processed_beds = []
        chrom_style = analyzer.get_chrom_style(extended_input)
        for i, bed in enumerate(bed_files):
            outb = os.path.join(temp_dir, f"bed_{i}.bed")
            analyzer.process_bed_file(bed, outb, chrom_style)
            processed_beds.append(outb)

        analyzer.run_bedtools_chain(extended_input, processed_beds, output_file)
        print(f"Overlap analysis complete. Output at {output_file}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
    return 0

def main():
    parser = argparse.ArgumentParser(prog="chromacs-overlap-expr",
                                     description="Overlap / expression helper. If no subcommand given, launches GUI.")
    sub = parser.add_subparsers(dest="cmd")

    ov = sub.add_parser("overlap", help="Run overlap analysis from CLI")
    ov.add_argument("--input", required=True, help="Primary ATAC-seq peaks BED")
    ov.add_argument("--bed", dest="beds", action="append", required=True, help="Additional BED file (can supply multiple)")
    ov.add_argument("--output", required=True, help="Output file")
    ov.add_argument("--window", type=int, default=200, help="Window size in bp")
    ov.set_defaults(func=cli_overlap)

    args, remaining = parser.parse_known_args()
    if args.cmd == "overlap":
        return_code = args.func(args)
        sys.exit(return_code)
    else:
        root = tk.Tk()
        app = ATACSeqAnalyzer(root)
        root.mainloop()

if __name__ == "__main__":
    main()
