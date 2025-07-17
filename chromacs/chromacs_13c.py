import tkinter as tk
import tkinter.font as tkFont
from tkinter import filedialog, messagebox
from tkinter import ttk
import subprocess
import os
import threading
import glob
import re
import pandas as pd
import sys
from pkg_resources import resource_filename
import gzip
import shutil
import webbrowser


# font setup
def load_custom_font(root):
    font_path = resource_filename("chromacs", "fonts/Roboto-Regular.ttf")

    if os.path.exists(font_path):
        try:
            roboto_font = tkFont.Font(root=root, family="Roboto", size=10)
            print(f"Successfully loaded custom font: {font_path}")
            return roboto_font

        except Exception as e:
            print(f"Error loading custom font: {e}")
            return None
    else:
        print(f"Custom font file not found at {font_path}. Using system default font.")
        return None


class ATACSeqPipeline:
    def __init__(self, root):
        self.root = root
        self.root.title("CHROMACS")
        self.root.geometry("950x800")

        # Top logo frame
        top_logo_frame = tk.Frame(self.root)
        top_logo_frame.pack(side="top", fill="x", pady=(5, 10), padx=10)

        # Left side: ChromAcS logo
        chromacs_logo_path = resource_filename("chromacs", "assets/ChromAcS.png")
        if os.path.exists(chromacs_logo_path):
            chromacs_logo_img = tk.PhotoImage(file=chromacs_logo_path).subsample(3, 3)
            chromacs_logo_label = tk.Label(top_logo_frame, image=chromacs_logo_img)
            chromacs_logo_label.image = chromacs_logo_img  # prevent garbage collection
            chromacs_logo_label.pack(side="left")

        # Right side: lab logo + lab name
        lab_frame = tk.Frame(top_logo_frame)
        lab_frame.pack(side="right")

        lab_logo_path = resource_filename("chromacs", "assets/lab_logo.png")
        if os.path.exists(lab_logo_path):
            lab_logo_img = tk.PhotoImage(file=lab_logo_path).subsample(3, 3)
            lab_logo_label = tk.Label(lab_frame, image=lab_logo_img)
            lab_logo_label.image = lab_logo_img  # prevent garbage collection
            lab_logo_label.pack(side="left", padx=(0, 5))

        def open_lab_website(event=None):
            webbrowser.open("https://www.epigen-bioinfolab.com/")

        lab_name_label = tk.Label(
            lab_frame,
            text="Epigen-BioinfoLab",
            font=("Roboto", 11, "bold"),
            fg="Green",
            cursor="hand2"
        )
        lab_name_label.pack(side="left")
        lab_name_label.bind("<Button-1>", open_lab_website)

        # Title bar icon (Linux .xbm format)
        icon_path = resource_filename("chromacs", "assets/ChromAcS.xbm")
        if os.path.exists(icon_path):
            self.root.iconbitmap(f"@{icon_path}")

        # Load custom font
        self.roboto_font = load_custom_font(self.root)

        if self.roboto_font:
            style = ttk.Style()
            style.configure('.', font=self.roboto_font)
        else:
            style = ttk.Style()
            style.configure('.', font=("TkDefaultFont", 10))

        # Dictionary to store parameters from each step
        self.params = {
            "step1": {},  # Raw data selection; then fastqc and multiqc
            "step2": {},  # Trim Galore (and sample names)
            "step3": {},  # Ref genome for alignment, indexing, and sorting
            "step4": {},  # Peak Calling (Genrich or MACS3)
            "step5": {}  # Peak Annotation, plus diffbind and noiseq on separate modules
        }

        # Creating a Notebook widget for the wizard steps
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True)

        # Creating frames for each step.
        self.step1_frame = ttk.Frame(self.notebook)
        self.step2_frame = ttk.Frame(self.notebook)
        self.step3_frame = ttk.Frame(self.notebook)
        self.step4_frame = ttk.Frame(self.notebook)
        self.step5_frame = ttk.Frame(self.notebook)

        # Adding frames (tabs) to the Notebook
        self.notebook.add(self.step1_frame, text="Step 1")
        self.notebook.add(self.step2_frame, text="Step 2")
        self.notebook.add(self.step3_frame, text="Step 3")
        self.notebook.add(self.step4_frame, text="Step 4")
        self.notebook.add(self.step5_frame, text="Step 5")

        # Disabling tabs for steps that are not yet completed

        self.notebook.tab(1, state="disabled")
        self.notebook.tab(2, state="disabled")
        self.notebook.tab(3, state="disabled")
        self.notebook.tab(4, state="disabled")

        # Building each step's UI
        self.setup_step1_ui()
        self.setup_step2_ui()
        self.setup_step3_ui()
        self.setup_step4_ui()
        self.setup_step5_ui()

        # Container to hold output + logos (bottom portion of the app)
        main_bottom_frame = tk.Frame(self.root)
        main_bottom_frame.pack(side="bottom", fill="both", expand=True)

        # Output frame inside bottom container
        output_frame = tk.Frame(main_bottom_frame)
        output_frame.pack(side="top", fill="both", expand=True, padx=10, pady=(10, 0))

        # Creating a Text widget for output (shows real-time run)
        self.output_text = tk.Text(output_frame, wrap=tk.WORD, height=20, width=100)
        self.output_text.pack(side="left", fill="both", expand=True)

        scrollbar = tk.Scrollbar(output_frame, command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")

        self.output_text.config(yscrollcommand=scrollbar.set)
        self.output_text.config(state=tk.DISABLED)

    def get_timestamp(self):
        import datetime
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def clear_frame(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()

    # =================== UI Setup for Each Step =================== #

    # Step 1: Selection of raw data ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def setup_step1_ui(self):
        # Output directory widgets
        tk.Label(self.step1_frame, text="Base Output Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=0,
                                                                                                            column=0,
                                                                                                            sticky="w",
                                                                                                            padx=10,
                                                                                                            pady=5)
        self.output_dir = tk.Entry(self.step1_frame, width=50)
        self.output_dir.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_output_dir, bg='grey').grid(row=0, column=2,
                                                                                                   padx=10, pady=5)

        instruction_text_2 = " [Assign the Base Output Directory where all the corresponding results will be saved. ]"
        tk.Label(self.step1_frame, text=instruction_text_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=0, columnspan=3, padx=10, pady=5, sticky="w")

        tk.Label(self.step1_frame, text="Raw Data Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=2,
                                                                                                         column=0,
                                                                                                         sticky="w",
                                                                                                         padx=10,
                                                                                                         pady=5)
        self.raw_data_dir = tk.Entry(self.step1_frame, width=50)
        self.raw_data_dir.grid(row=2, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_raw_data, bg='grey').grid(row=2, column=2,
                                                                                                 padx=10, pady=5)

        instruction_text_1 = "[ If 'Select specific samples' is unchecked, ALL FASTQ files in the RAW DATA DIRECTORY will be processed. ]"
        tk.Label(self.step1_frame, text=instruction_text_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=0, columnspan=3, padx=10, pady=5, sticky="w")

        # Checkbox for selective sample selection
        self.select_samples_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            self.step1_frame,
            text="Select specific samples",
            variable=self.select_samples_var,
            command=self.toggle_sample_selection
        ).grid(row=4, column=0, columnspan=3, pady=5)

        # Listbox for sample selection
        self.sample_listbox = tk.Listbox(self.step1_frame, selectmode=tk.MULTIPLE, height=6, width=70)
        self.sample_listbox.grid(row=5, column=0, columnspan=3, padx=10, pady=5)
        self.sample_listbox.grid_remove()

        # Scrollbar for listbox
        scrollbar = tk.Scrollbar(self.step1_frame, orient="vertical")
        scrollbar.config(command=self.sample_listbox.yview)
        scrollbar.grid(row=5, column=3, sticky="ns")
        self.sample_listbox.config(yscrollcommand=scrollbar.set)

        tk.Label(self.step1_frame, text="Number of Threads (Default- 8):", font=(self.roboto_font, 10, 'bold')).grid(
            row=6, column=0, sticky="w", padx=10, pady=5)
        self.threads_entry = tk.Entry(self.step1_frame, width=10)
        self.threads_entry.grid(row=6, column=1, padx=10, pady=5)
        self.threads_entry.insert(0, "8")  # Default value

        # Save button
        tk.Button(self.step1_frame, text="Save & Next", command=self.save_step1_next, bg="yellow green").grid(
            row=7, column=0, columnspan=3, pady=10)

        # helper functions of step1_ui

    def browse_raw_data(self):
        directory = filedialog.askdirectory()
        if directory:
            self.raw_data_dir.delete(0, tk.END)
            self.raw_data_dir.insert(0, directory)
            # Auto-populate sample list if checkbox is checked
            if self.select_samples_var.get():
                self.populate_sample_list()

    def browse_output_dir(self):
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir.delete(0, tk.END)
            self.output_dir.insert(0, directory)

    def toggle_sample_selection(self):
        if self.select_samples_var.get():
            self.sample_listbox.grid()
            self.populate_sample_list()
        else:
            self.sample_listbox.grid_remove()

    def populate_sample_list(self):
        self.sample_listbox.delete(0, tk.END)
        directory = self.raw_data_dir.get()
        if directory:
            # Find all FASTQ files with both extensions (*fastq.gz and *fq.gz)
            fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz")) + glob.glob(
                os.path.join(directory, "*.fq.gz"))
            for f in sorted(fastq_files):
                self.sample_listbox.insert(tk.END, os.path.basename(f))

    def save_step1_next(self):
        self.params["step1"]["raw_data_dir"] = self.raw_data_dir.get().strip()
        self.params["step1"]["base_output_dir"] = self.output_dir.get().strip()

        # Store selected samples WITH FULL PATHS if selection is enabled
        if self.select_samples_var.get():
            raw_data_dir = self.params["step1"]["raw_data_dir"]
            # Get selected filenames from listbox and convert to full paths
            selected_files = [
                os.path.join(raw_data_dir, self.sample_listbox.get(i))
                for i in self.sample_listbox.curselection()
            ]
            self.params["step1"]["selected_samples"] = selected_files
        else:
            self.params["step1"]["selected_samples"] = None  # Process all

        if not self.params["step1"]["raw_data_dir"] or not self.params["step1"]["base_output_dir"]:
            messagebox.showerror("Error", "Please fill in both directories.")
            return

        # Extract base sample names for Step 2 display
        sample_names = set()
        raw_data_dir = self.params["step1"]["raw_data_dir"]
        selected_files = self.params["step1"].get("selected_samples", [])

        # If no specific selection, find all FASTQ files (move this inside the above condition where else = None)
        if not selected_files:
            selected_files = glob.glob(os.path.join(raw_data_dir, "*.fastq.gz")) + \
                             glob.glob(os.path.join(raw_data_dir, "*.fq.gz"))

        # Extract base sample names from full paths
        for filepath in selected_files:
            base = os.path.basename(filepath)

            # Match: _R1_001, _1.fastq.gz, etc.
            match = re.match(r"^(.*?)(?:_R?[12](?:_001)?|_[12])\.f(?:ast)?q\.gz$", base)

            if match:
                sample_name = match.group(1)
            else:
                # Fallback: strip known suffixes
                sample_name = re.sub(r'(_R?[12](?:_001)?|_[12])\.f(?:ast)?q\.gz$', '', base)

            sample_names.add(sample_name)

        self.params["step1"]["auto_sample_names"] = sorted(sample_names)

        # update thread as set by the user
        threads = self.threads_entry.get().strip()
        try:
            threads = int(threads)
            if threads < 1:
                raise ValueError("Thread count must be >= 1")
        except ValueError:
            threads = 8  # fallback
            self.update_output_gui("Invalid thread count. Defaulting to 8.\n")

        self.params["threads"] = threads

        # Update Step 2 UI
        self.sample_listbox_step2.delete(0, tk.END)
        for name in self.params["step1"]["auto_sample_names"]:
            self.sample_listbox_step2.insert(tk.END, name)

        self.notebook.tab(1, state="normal")
        self.notebook.select(1)

    # Step 2: Trim Galore and Sample Names +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def setup_step2_ui(self):
        tk.Label(self.step2_frame, text="Auto-Detected Sample Names:", font=(self.roboto_font, 10, 'bold')).grid(
            row=0, column=0, sticky="w", padx=10, pady=5)

        # Listbox to display samples (read-only)
        self.sample_listbox_step2 = tk.Listbox(self.step2_frame, selectmode=tk.SINGLE, height=6, width=70)
        self.sample_listbox_step2.grid(row=0, column=1, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Samples auto-extracted from filenames.\n"
            "Format: SRRXXXXX (ignores _1/_2 suffixes).\n\n"
        )

        tk.Label(self.step2_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=1, columnspan=3, padx=10, pady=5, sticky="w")

        self.skip_trimming_var = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step2_frame, text="Skip Trimming (Use Raw Data for Alignment)",
                       variable=self.skip_trimming_var, font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0,
                                                                                                  columnspan=2, padx=10,
                                                                                                  pady=5)

        # Instructional text
        instruction_3 = (
            "Checking this will skip the tirmming done by Trim Galore.\n"
            "Raw data will be used directly for alignment.\n"
            "Please check the FastQC and MultiQC reports if you are unsure of this step. \n"
        )

        tk.Label(self.step2_frame, text=instruction_3, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=5, column=1, columnspan=3, padx=10, pady=5, sticky="w")

        tk.Button(self.step2_frame, text="Save & Next", command=self.save_step2_next, bg="yellow green").grid(row=7,
                                                                                                              column=1,
                                                                                                              pady=5)
        tk.Button(self.step2_frame, text="Back", command=lambda: self.notebook.select(0), bg='salmon').grid(row=7,
                                                                                                            column=0,
                                                                                                            pady=5)

        # helper functions of step2_ui

    def save_step2_next(self):
        # Get auto-detected names from Step 1
        self.params["step2"]["sample_names"] = self.params["step1"]["auto_sample_names"]
        self.params["step2"]["skip_trimming"] = self.skip_trimming_var.get()

        if not self.params["step2"]["sample_names"]:
            messagebox.showerror("Error", "No samples detected. Check Step 1.")
            return

        self.notebook.tab(2, state="normal")
        self.notebook.select(2)

    # Step 3: Alignment Settings +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def setup_step3_ui(self):
        # Genome mapping with valid MACS3 parameters
        self.genome_map = {
            # === Vertebrates ===
            "Human (GRCh38)": {
                "bowtie_ref": "GRCh38",
                "macs_size": "hs",
                "ensembl_species": "homo_sapiens",
                "ensembl_cap": "Homo_sapiens",
                "ensembl_assembly": "GRCh38"
            },
            "Mouse (GRCm39)": {
                "bowtie_ref": "GRCm39",
                "macs_size": "mm",
                "ensembl_species": "mus_musculus",
                "ensembl_cap": "Mus_musculus",
                "ensembl_assembly": "GRCm39"
            },
            "Rat (mRatBN7.2)": {
                "bowtie_ref": "mRatBN7.2",
                "macs_size": "2.9e9",
                "ensembl_species": "rattus_norvegicus",
                "ensembl_cap": "Rattus_norvegicus",
                "ensembl_assembly": "mRatBN7.2"
            },
            "Cow (ARS-UCD1.3)": {
                "bowtie_ref": "ARS-UCD1.3",
                "macs_size": "2.7e9",
                "ensembl_species": "bos_taurus",
                "ensembl_cap": "Bos_taurus",
                "ensembl_assembly": "ARS-UCD1.3"
            },
            "Pig (Sscrofa11.1)": {
                "bowtie_ref": "Sscrofa11.1",
                "macs_size": "2.5e9",
                "ensembl_species": "sus_scrofa",
                "ensembl_cap": "Sus_scrofa",
                "ensembl_assembly": "Sscrofa11.1"
            },
            "Chicken (GRCg7b)": {
                "bowtie_ref": "GRCg7b",
                "macs_size": "1.2e9",
                "ensembl_species": "gallus_gallus",
                "ensembl_cap": "Gallus_gallus",
                "ensembl_assembly": "bGalGal1.mat.broiler.GRCg7b"
            },
            "Chimpanzee (Pan_tro_3.0)": {
                "bowtie_ref": "Pan_tro_3.0",
                "macs_size": "3.3e9",
                "ensembl_species": "pan_troglodytes",
                "ensembl_cap": "Pan_troglodytes",
                "ensembl_assembly": "Pan_tro_3.0"
            },

            "Dog (ROS_Cfam_1.0)": {
                "bowtie_ref": "ROS_Cfam_1.0",
                "macs_size": "2.4e9",
                "ensembl_species": "canis_lupus_familiaris",
                "ensembl_cap": "Canis_lupus_familiaris",
                "ensembl_assembly": "ROS_Cfam_1.0"
            },
            "Goat (ARS1)": {
                "bowtie_ref": "ARS1",
                "macs_size": "2.9e9",
                "ensembl_species": "capra_hircus",
                "ensembl_cap": "Capra_hircus",
                "ensembl_assembly": "ARS1"
            },
            "Rabbit (OryCun2.0)": {
                "bowtie_ref": "OryCun2.0",
                "macs_size": "2.7e9",
                "ensembl_species": "oryctolagus_cuniculus",
                "ensembl_cap": "Oryctolagus_cuniculus",
                "ensembl_assembly": "OryCun2.0"
            },
            "Gorilla (gorGor4)": {
                "bowtie_ref": "gorGor4",
                "macs_size": "3.5e9",
                "ensembl_species": "gorilla_gorilla",
                "ensembl_cap": "Gorilla_gorilla",
                "ensembl_assembly": "gorGor4"
            },
            "Rhesus Macaque (Mmul_10)": {
                "bowtie_ref": "Mmul_10",
                "macs_size": "2.9e9",
                "ensembl_species": "macaca_mulatta",
                "ensembl_cap": "Macaca_mulatta",
                "ensembl_assembly": "Mmul_10"
            },

            # === Fish & Amphibian ===
            "Zebrafish (GRCz11)": {
                "bowtie_ref": "GRCz11",
                "macs_size": "1.4e9",
                "ensembl_species": "danio_rerio",
                "ensembl_cap": "Danio_rerio",
                "ensembl_assembly": "GRCz11"
            },
            "Frog (UCB_Xtro_10.0)": {
                "bowtie_ref": "UCB_Xtro_10.0",
                "macs_size": "1.7e9",
                "ensembl_species": "xenopus_tropicalis",
                "ensembl_cap": "Xenopus_tropicalis",
                "ensembl_assembly": "UCB_Xtro_10.0"
            },
            "Indian Cobra (Nana_v5)": {
                "bowtie_ref": "Nana_v5",
                "macs_size": "1.8e9",
                "ensembl_species": "naja_naja",
                "ensembl_cap": "Naja_naja",
                "ensembl_assembly": "Nana_v5"
            },
            "Atlantic Salmon (Ssal_v3.1)": {
                "bowtie_ref": "Ssal_v3.1",
                "macs_size": "3.0e9",
                "ensembl_species": "salmo_salar",
                "ensembl_cap": "Salmo_salar",
                "ensembl_assembly": "Ssal_v3.1"
            },
            # === Invertebrates ===
            "C. elegans (WBcel235)": {
                "bowtie_ref": "WBcel235",
                "macs_size": "1.0e8",
                "ensembl_species": "caenorhabditis_elegans",
                "ensembl_cap": "Caenorhabditis_elegans",
                "ensembl_assembly": "WBcel235"
            },
            "Drosophila (BDGP6)": {
                "bowtie_ref": "BDGP6",
                "macs_size": "1.4e8",
                "ensembl_species": "drosophila_melanogaster",
                "ensembl_cap": "Drosophila_melanogaster",
                "ensembl_assembly": "BDGP6"
            },
        }

        # Genome selection
        tk.Label(self.step3_frame, text="Select Organism/Genome:", font=(self.roboto_font, 10, 'bold')).grid(row=0,
                                                                                                             column=0,
                                                                                                             sticky="w",
                                                                                                             padx=10,
                                                                                                             pady=5)
        self.genome_var = tk.StringVar()
        self.genome_menu = ttk.Combobox(self.step3_frame, textvariable=self.genome_var,
                                        values=list(self.genome_map.keys()), state="readonly")
        self.genome_menu.grid(row=0, column=1, padx=10, pady=5)

        # Reference status indicator
        self.ref_status = tk.Label(self.step3_frame, text="", fg="salmon")
        self.ref_status.grid(row=1, column=1, sticky="w", padx=10, pady=5)

        # Buttons
        tk.Button(self.step3_frame, text="Check Reference", command=self.check_reference, bg="gray").grid(row=2,
                                                                                                          column=1,
                                                                                                          pady=5)
        tk.Button(self.step3_frame, text="Save & Next", command=self.save_step3_next, bg="yellow green").grid(row=3,
                                                                                                              column=1,
                                                                                                              pady=5)
        tk.Button(self.step3_frame, text="Back", command=lambda: self.notebook.select(1), bg="salmon").grid(row=3,
                                                                                                            column=0,
                                                                                                            pady=5)

    def check_reference(self):
        selected = self.genome_var.get()
        if not selected:
            messagebox.showerror("Error", "Please select a genome reference")
            return False

        genome_info = self.genome_map[selected]

        # Use the script's directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        ref_dir = os.path.join(script_dir, "ref_genome", genome_info["bowtie_ref"])

        bt2_base = os.path.join(ref_dir, genome_info["bowtie_ref"])
        fa_file = os.path.join(ref_dir, f"{genome_info['bowtie_ref']}.fa")
        fa_gz = f"{fa_file}.gz"

        # Check index files
        index_files = [f"{bt2_base}.{ext}" for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]
        index_exists = all(os.path.exists(f) for f in index_files)

        if index_exists:
            self.ref_status.config(text="Reference index found!", fg="yellow green")
            return True
        else:
            status_messages = []
            if os.path.exists(fa_file):
                status_messages.append("FASTA file exists")
            elif os.path.exists(fa_gz):
                status_messages.append("Compressed FASTA exists")
            else:
                status_messages.append("No reference files found")

            status_text = "Reference status: " + ", ".join(status_messages) + ". Will be processed during pipeline run."
            self.ref_status.config(text=status_text, fg="orange")

            # Save everything in params
            self.params["step3"].update({
                "genome_version": genome_info["bowtie_ref"],
                "ref_dir": ref_dir,
                "bt2_base": bt2_base,
                "ensembl_species": genome_info["ensembl_species"],
                "ensembl_cap": genome_info["ensembl_cap"],
                "ensembl_assembly": genome_info["ensembl_assembly"]
            })
            return False

    def save_step3_next(self):
        if not self.genome_var.get():
            messagebox.showerror("Error", "Select a genome reference.")
            return

        selected = self.genome_var.get()
        genome_info = self.genome_map[selected]

        script_dir = os.path.dirname(os.path.abspath(__file__))
        ref_dir = os.path.join(script_dir, "ref_genome", genome_info["bowtie_ref"])
        bt2_base = os.path.join(ref_dir, genome_info["bowtie_ref"])

        self.params["step3"].update({
            "genome_version": genome_info["bowtie_ref"],
            "macs3_genome_size": genome_info["macs_size"],
            "ensembl_species": genome_info["ensembl_species"],
            "ensembl_cap": genome_info["ensembl_cap"],
            "ensembl_assembly": genome_info["ensembl_assembly"],
            "ref_dir": ref_dir,
            "bt2_base": bt2_base
        })
        # Enable next step
        self.notebook.tab(3, state="normal")
        self.notebook.select(3)

    # Step 4: Peak Calling selection +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    def setup_step4_ui(self):
        tk.Label(self.step4_frame, text="Peak Calling", font=("Arial", 16, "bold")).grid(row=0, column=0,
                                                                                         columnspan=2, pady=10)
        tk.Label(self.step4_frame, text="Select Peak Caller:", font=(self.roboto_font, 10, 'bold')).grid(row=1,
                                                                                                         column=0,
                                                                                                         sticky="w",
                                                                                                         padx=10,
                                                                                                         pady=5)
        self.peak_caller_choice = tk.StringVar(value="Genrich")
        tk.Radiobutton(self.step4_frame, text="Genrich", variable=self.peak_caller_choice, value="Genrich").grid(row=2,
                                                                                                                 column=0,
                                                                                                                 sticky="w",
                                                                                                                 padx=10)
        tk.Radiobutton(self.step4_frame, text="MACS3", variable=self.peak_caller_choice, value="MACS3").grid(row=2,
                                                                                                             column=1,
                                                                                                             sticky="w",
                                                                                                             padx=10)
        tk.Button(self.step4_frame, text="Save & Next", command=self.save_step4_next, bg="yellow green").grid(row=3,
                                                                                                              column=2,
                                                                                                              columnspan=3,
                                                                                                              pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.notebook.select(2), bg='salmon').grid(row=3,
                                                                                                            column=0,
                                                                                                            pady=10)

    def save_step4_next(self):
        self.params["step4"]["peak_caller"] = self.peak_caller_choice.get().strip()
        # Depending on the choice, user can navigate to the corresponding sub-UI.
        if self.peak_caller_choice.get() == "Genrich":
            self.setup_step4a_genrich_ui()
        elif self.peak_caller_choice.get() == "MACS3":
            self.setup_step4b_macs3_ui()

    def setup_step4a_genrich_ui(self, samples=None):
        # Genrich UI (for Step 4a)
        self.clear_frame(self.step4_frame)
        tk.Label(self.step4_frame, text="Genrich Peak Calling", font=("Arial", 16, "bold")).grid(
            row=0, column=0, columnspan=3, pady=10)
        tk.Label(self.step4_frame, text="Select Peak Type:", font=(self.roboto_font, 10, 'bold')).grid(row=1, column=0,
                                                                                                       sticky="w",
                                                                                                       padx=10, pady=5)
        self.peak_type_choice = tk.StringVar(value="merged")
        self.peak_type_choice.trace("w", lambda *args: self.toggle_merged_output_field())
        tk.Radiobutton(self.step4_frame, text="Merged", variable=self.peak_type_choice, value="merged").grid(
            row=1, column=1, sticky="w", padx=10)
        tk.Radiobutton(self.step4_frame, text="Unmerged", variable=self.peak_type_choice, value="unmerged").grid(
            row=1, column=2, sticky="w", padx=10)

        # Instructional text
        instruction_1 = (
            "If UNMERGED is checked, SEPARATE PEAK FILES are generated for each sample from the pool of multiple samples\n"
            "If MERGED is checked, ONE MERGED PEAK FILE is generated from the pool of multiple samples [suitable for replicates]"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=1, column=3, columnspan=4, padx=10, pady=5, sticky="w")

        # Uses sample names from Step 2.
        step2_samples = self.params["step2"].get("sample_names", [])
        if not step2_samples:
            tk.Label(self.step4_frame, text="No sample names found. Please set them in Step 2.").grid(
                row=2, column=0, columnspan=3, padx=10, pady=5)
            return

        # Test samples
        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0,
                                                                                                   sticky="w", padx=10,
                                                                                                   pady=5)
        self.treated_sample_var = tk.StringVar(value="Select a test sample")
        treated_dropdown = tk.OptionMenu(self.step4_frame, self.treated_sample_var, *step2_samples)
        treated_dropdown.grid(row=3, column=1, padx=10, pady=5)
        self.treated_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.treated_listbox.grid(row=3, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Test Sample", command=self.add_treated_sample, bg='gray').grid(
            row=3, column=3, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Select one by one: click on a sample and then click on \'Add Test Sample\'\n"
            "Repeat this for the next sample \n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Control samples
        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0,
                                                                                                      sticky="w",
                                                                                                      padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=4, column=1, padx=10, pady=5)
        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=4, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(
            row=4, column=3, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put untreated condition as control here\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=4, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Blacklist files
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=5,
                                                                                                             column=0,
                                                                                                             sticky="w",
                                                                                                             padx=10,
                                                                                                             pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=5, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=5,
                                                                                                           column=2,
                                                                                                           padx=10,
                                                                                                           pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=5,
                                                                                                                 column=3,
                                                                                                                 padx=10,
                                                                                                                 pady=5)

        # Instructional text
        instruction_1 = (
            "Note: You can download blacklist from our Github page\n"
            "You can skip this part if you are unable to obtain a blacklist for your target organism\n"
            "Please ensure that you are using blacklists from the appropriate assembly versions\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=5, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Exclude chromosomes
        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=6,
                                                                                                          column=0,
                                                                                                          sticky="w",
                                                                                                          padx=10,
                                                                                                          pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (mitochondrial chromosomes)",
                       variable=self.exclude_chr_m).grid(row=6, column=1,
                                                         padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=6, column=2,
                                                                                             padx=10, pady=5)
        # chromosome exclusion entry
        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=6, column=3, sticky="w", padx=10,
                                                                                 pady=5)
        self.custom_chr_entry = tk.Entry(self.step4_frame, width=30)
        self.custom_chr_entry.grid(row=6, column=4, padx=10, pady=5)

        # expand_cut_sites
        tk.Label(self.step4_frame, text='Expand cut-sites to _ bp:', font=(self.roboto_font, 10, 'bold')).grid(row=7,
                                                                                                               column=0,
                                                                                                               padx=10,
                                                                                                               pady=5,
                                                                                                               sticky='w')
        self.expand_cut_sites_var = tk.StringVar(value="100")
        self.expand_cut_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.expand_cut_sites_var,
            values=["100", "120", "150"]
        )
        self.expand_cut_combobox.grid(row=7, column=1, padx=10, pady=5)

        # Instructional text
        instruction_2 = (
            "For custom chromosome entry ensure to enter the correct chromosome naming format, like for Homo sapiens: 22\n"
            "For multiple entries they must be SINGLE WHITE-SPACED; like, 1 2 X 22\n"
            "Do NOT use prefix like 'chr'"
        )

        tk.Label(self.step4_frame, text=instruction_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=7, column=4, columnspan=6, padx=10, pady=5, sticky="w")

        # max_q_value
        tk.Label(self.step4_frame, text='Assign Max q-value (0 to 1):', font=(self.roboto_font, 10, 'bold')).grid(row=8,
                                                                                                                  column=0,
                                                                                                                  padx=10,
                                                                                                                  pady=5)
        self.max_q_value_var = tk.StringVar(value="0.05")
        self.max_q_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.max_q_value_var,
            values=["0.05", "1"]
        )
        self.max_q_combobox.grid(row=8, column=1, padx=10, pady=5)

        # Merged output name
        self.merged_output_label = tk.Label(self.step4_frame, text="Merged Output Name:",
                                            font=(self.roboto_font, 10, 'bold'))
        self.merged_output_label.grid(row=9, column=0, sticky="w", padx=10, pady=5)
        self.merged_output_name_var = tk.StringVar()
        self.merged_output_entry = tk.Entry(self.step4_frame, textvariable=self.merged_output_name_var, width=30)
        self.merged_output_entry.grid(row=9, column=1, padx=10, pady=5)

        # Save button (assign to self.save_button)
        self.save_button = tk.Button(
            self.step4_frame,
            text="Save and Next",
            command=self.save_step4a_settings,
            bg="yellow green"
        )
        self.save_button.grid(row=12, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='salmon').grid(
            row=12,
            column=0,
            pady=10)

    def save_step4a_settings(self):
        self.params["step4"]["peak_type"] = self.peak_type_choice.get().strip()
        self.params["step4"]["treated_samples"] = [self.treated_listbox.get(i) for i in
                                                   self.treated_listbox.curselection()]
        self.params["step4"]["control_samples"] = [self.control_listbox.get(i) for i in
                                                   self.control_listbox.curselection()]
        self.params["step4"]["blacklist_files"] = [self.blacklist_listbox.get(i) for i in
                                                   range(self.blacklist_listbox.size())]

        # Process exclude chromosomes:
        exclude_chr_list = []
        if self.exclude_chr_m.get():
            exclude_chr_list.append("MT")
        if self.exclude_chr_y.get():
            exclude_chr_list.append("Y")
        # custom input as space-separated tokens
        custom_chr = self.custom_chr_entry.get().strip()
        if custom_chr:
            # split on whitespace
            custom_list = custom_chr.split()
            exclude_chr_list.extend(custom_list)
        self.params["step4"]["exclude_chr"] = exclude_chr_list

        # Process numeric parameters:
        try:
            expand_cut = int(self.expand_cut_sites_var.get().strip())
            max_q = float(self.max_q_value_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric input for expand cut-sites or Q-value!")
            return

        self.params["step4"]["expand_cut_sites"] = expand_cut
        self.params["step4"]["max_q_value"] = max_q
        self.params["step4"]["merged_output_name"] = self.merged_output_name_var.get().strip()

        self.notebook.tab(4, state="normal")
        self.notebook.select(4)
        messagebox.showinfo("Step 4a Saved", "Genrich settings saved.")

    def toggle_merged_output_field(self):
        if self.peak_type_choice.get() == "merged":
            # Show merged output widgets at fixed rows
            self.merged_output_label.grid(row=9, column=0, sticky="w", padx=10, pady=5)
            self.merged_output_entry.grid(row=9, column=1, padx=10, pady=5)
            # Move the Save button down to row 10
            self.save_button.grid(row=10, column=0, columnspan=3, pady=10)
        else:
            # Hide merged output widgets
            self.merged_output_label.grid_remove()
            self.merged_output_entry.grid_remove()
            # Move the Save button up to row 9
            self.save_button.grid(row=12, column=0, columnspan=3, pady=10)

    def add_treated_sample(self):
        sample = self.treated_sample_var.get()
        if sample and sample != "Select a test sample" and sample not in self.treated_listbox.get(0, tk.END):
            self.treated_listbox.insert(tk.END, sample)
            self.treated_listbox.select_set(tk.END)
        else:
            print("Test sample not added.")

    def add_control_sample(self):
        sample = self.control_sample_var.get()
        if sample and sample != "Select a control sample" and sample not in self.control_listbox.get(0, tk.END):
            self.control_listbox.insert(tk.END, sample)
            self.control_listbox.select_set(tk.END)
        else:
            print("Control sample not added.")

    def browse_blacklist_files(self):
        files = filedialog.askopenfilenames(
            title="Select Blacklist Files",
            filetypes=[("BED files", "*.bed"), ("All Files", "*.*")]
        )
        for file in files:
            if file not in self.blacklist_listbox.get(0, tk.END):
                self.blacklist_listbox.insert(tk.END, file)

    def remove_blacklist_files(self):
        selected_indices = self.blacklist_listbox.curselection()
        for idx in reversed(selected_indices):
            self.blacklist_listbox.delete(idx)

    def navigate_back_to_step4_ui(self):
        self.clear_frame(self.step4_frame)
        self.setup_step4_ui()

    def setup_step4b_macs3_ui(self, samples=None):
        # MACS3 UI (for Step 4b).
        self.clear_frame(self.step4_frame)
        tk.Label(self.step4_frame, text="MACS3 Peak Calling", font=("Arial", 16, "bold")).grid(
            row=0, column=0, columnspan=3, pady=10)

        # uses sample names from step2
        if samples is None:
            step2_samples = self.params["step2"].get("sample_names", [])
        else:
            step2_samples = samples
        if not step2_samples:
            tk.Label(self.step4_frame, text="Enter Sample Names (separated by spaces):").grid(
                row=1, column=0, sticky="w", padx=10, pady=5)
            self.sample_input = tk.Entry(self.step4_frame, width=50)
            self.sample_input.grid(row=1, column=1, columnspan=2, padx=10, pady=5)
            tk.Button(self.step4_frame, text="Submit Samples", command=self.populate_samples).grid(
                row=1, column=3, padx=10, pady=5)
            return

        # Test samples dropdown + Listbox
        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=2, column=0,
                                                                                                   sticky="w", padx=10,
                                                                                                   pady=5)
        self.treated_sample_var = tk.StringVar(value="Select a test sample")
        treated_dropdown = tk.OptionMenu(self.step4_frame, self.treated_sample_var, *step2_samples)

        treated_dropdown.grid(row=2, column=1, padx=10, pady=5)
        self.treated_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.treated_listbox.grid(row=2, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Test Sample", command=self.add_treated_sample, bg='gray').grid(
            row=2, column=3, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Select one by one: click on a sample and then click on \'Add Test Sample\'\n"
            "Repeat this for the next sample; you can ADD MULTIPLE SAMPLES at a time \n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=2, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Control samples dropdown + Listbox
        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0,
                                                                                                      sticky="w",
                                                                                                      padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=3, column=1, padx=10, pady=5)

        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=3, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(
            row=3, column=3, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put untreated condition as control here\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Blacklist file selection
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=4,
                                                                                                             column=0,
                                                                                                             sticky="w",
                                                                                                             padx=10,
                                                                                                             pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=4, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=4,
                                                                                                           column=2,
                                                                                                           padx=10,
                                                                                                           pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=4,
                                                                                                                 column=3,
                                                                                                                 padx=10,
                                                                                                                 pady=5)

        # Instructional text
        instruction_1 = (
            "Note: You can download blacklist from our Github page\n"
            "You can skip this part if you are unable to obtain a blacklist for your target organism\n"
            "Please ensure that you are using blacklists from the appropriate assembly versions\n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=4, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Options for excluding chromosomes
        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=5,
                                                                                                          column=0,
                                                                                                          sticky="w",
                                                                                                          padx=10,
                                                                                                          pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (Mitochondrial chromosomes)",
                       variable=self.exclude_chr_m).grid(row=5, column=1,
                                                         padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=5, column=2,
                                                                                             padx=10, pady=5)
        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=5, column=3, sticky="w", padx=10,
                                                                                 pady=5)
        self.custom_chr_entry = tk.Entry(self.step4_frame, width=30)
        self.custom_chr_entry.grid(row=5, column=4, padx=10, pady=5)

        # Q-value selection
        tk.Label(self.step4_frame, text='Assign Max q-value (0 to 1):', font=(self.roboto_font, 10, 'bold')).grid(row=6,
                                                                                                                  column=0,
                                                                                                                  padx=10,
                                                                                                                  pady=5,
                                                                                                                  sticky="w")
        self.max_q_value_var = tk.StringVar(value="0.05")
        self.max_q_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.max_q_value_var,
            values=["0.05", "0.01"]
        )
        self.max_q_combobox.grid(row=6, column=1, padx=10, pady=5)

        # Instructional text
        instruction_2 = (
            "For custom chromosome entry ensure to enter the correct chromosome naming format, like for Homo sapiens: 22\n"
            "For multiple entries they must be SINGLE WHITE-SPACED; like, 1 2 X 22\n"
            "Do NOT use prefix like 'chr'"
        )

        tk.Label(self.step4_frame, text=instruction_2, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=6, column=4, columnspan=6, padx=10, pady=5, sticky="w")

        # With read-only display
        tk.Label(self.step4_frame, text="Genome Size (auto-filled):").grid(row=7, column=0, sticky="w", padx=10, pady=5)
        self.macs3_genome_entry = tk.Entry(self.step4_frame, width=30, state='readonly')
        self.macs3_genome_entry.grid(row=7, column=1, padx=10, pady=5)

        # Auto-populate from Step 3
        genome_size = self.params["step3"].get("macs3_genome_size", "hs")
        self.macs3_genome_entry.config(state='normal')
        self.macs3_genome_entry.insert(0, genome_size)
        self.macs3_genome_entry.config(state='readonly')
        tk.Button(self.step4_frame, text="Save and Next", command=self.save_step4b_settings, bg="yellow green").grid(
            row=9, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='salmon').grid(
            row=9,
            column=0,
            pady=10)

    def populate_samples(self):
        sample_input = self.sample_input.get()
        samples = sample_input.split()
        self.setup_step4b_macs3_ui(samples)

    def save_step4b_settings(self):
        self.params["step4"]["treated_samples"] = [self.treated_listbox.get(i) for i in
                                                   self.treated_listbox.curselection()]
        self.params["step4"]["control_samples"] = [self.control_listbox.get(i) for i in
                                                   self.control_listbox.curselection()]
        self.params["step4"]["blacklist_files"] = [self.blacklist_listbox.get(i) for i in
                                                   range(self.blacklist_listbox.size())]
        exclude_chr_list = []
        if self.exclude_chr_m.get():
            exclude_chr_list.append("MT")
        if self.exclude_chr_y.get():
            exclude_chr_list.append("Y")
        custom_chr = self.custom_chr_entry.get().strip()
        if custom_chr:
            custom_list = custom_chr.split()  # Split on whitespace
            exclude_chr_list.extend(custom_list)
        self.params["step4"]["exclude_chr"] = exclude_chr_list

        try:
            max_q = float(self.max_q_value_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Invalid Q-value!")
            return

        self.params["step4"]["max_q_value"] = max_q

        self.params["step4"]["genome"] = self.params["step3"]["macs3_genome_size"]

        messagebox.showinfo("Step 4b Saved", "MACS3 settings saved.")

        self.notebook.tab(4, state="normal")
        self.notebook.select(4)

    def setup_step5_ui(self):
        frame = self.step5_frame

        # Configure grid for responsiveness
        for i in range(3):  # 3 columns: left spacer, center, right spacer
            frame.grid_columnconfigure(i, weight=1)
        frame.grid_rowconfigure(0, weight=1)

        # Instruction label (centered)
        instruction_1 = (
            "Click 'Run Pipeline' to generate peak files with annotation.\n"
            "The first run may take longer due to reference genome setup and indexing.\n\n"
            "- FastQC, MultiQC, Trimmed data, Aligned BAM files, and Coverage files will be stored under the base output directory.\n"
            "- Peak files will be saved in the 'peak_files' folder.\n"
            "- Annotated peaks will be in 'Annotated_Peaks' inside 'peak_files'."
        )
        tk.Label(
            frame, text=instruction_1, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 10, 'bold')
        ).grid(row=0, column=1, padx=20, pady=10, sticky="ew")

        # Button row: Back and Run Pipeline (centered)
        button_frame = tk.Frame(frame)
        button_frame.grid(row=1, column=1, pady=10)

        tk.Button(button_frame, text="Back", command=lambda: self.notebook.select(3), bg='salmon'
                  ).pack(side=tk.LEFT, padx=20)
        tk.Button(button_frame, text="Run Pipeline", command=self.run_pipeline, bg="yellow green"
                  ).pack(side=tk.LEFT, padx=20)

        # DiffBind info (centered and responsive)
        diffbind_info = (
            "After peak calling is complete, you can run DiffBind to perform differential binding analysis "
            "based on experimental conditions."
        )
        tk.Label(
            frame, text=diffbind_info, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 9, 'italic')
        ).grid(row=2, column=1, padx=20, pady=(20, 5), sticky="ew")

        # DiffBind button (centered)
        diffbind_btn_frame = tk.Frame(frame)
        diffbind_btn_frame.grid(row=3, column=1, pady=10)
        self.diffbind_btn = tk.Button(
            diffbind_btn_frame, text="Run DiffBind Analysis", bg='orange',
            command=self.launch_diffbind_config, state=tk.NORMAL
        )
        self.diffbind_btn.pack()

        # NOISeq info label
        noisq_info = (
            "If you do not have biological replicates, you can use NOISeq instead of DiffBind.\n"
            "NOISeq performs a simulation-based differential analysis on a peak count matrix."
        )
        tk.Label(
            frame, text=noisq_info, wraplength=800, justify=tk.LEFT,
            anchor="center", font=(self.roboto_font, 9, 'italic')
        ).grid(row=4, column=1, padx=20, pady=(10, 5), sticky="ew")

        # NOISeq button (centered below)
        noisq_btn_frame = tk.Frame(frame)
        noisq_btn_frame.grid(row=5, column=1, pady=(5, 20))
        self.noisq_btn = tk.Button(
            noisq_btn_frame, text="Run NOISeq Analysis (No Replicates)", bg='light blue',
            command=self.launch_noisq_config
        )
        self.noisq_btn.pack()

        # Let rows expand if window height increases
        for i in range(5):
            frame.grid_rowconfigure(i, weight=1)

    # =================== Final Run Pipeline =================== #

    def run_pipeline(self):
        base_out = self.params["step1"].get("base_output_dir")
        if not base_out:
            messagebox.showerror("Error", "Base output directory not specified in Step 1!")
            return

        fastqc_raw = os.path.join(base_out, "fastqc_raw")
        multiqc_raw = os.path.join(base_out, "multiqc_raw")
        trimmed_data = os.path.join(base_out, "trimmed_data")
        bam_output = os.path.join(base_out, "bam_output")
        fastqc_trimmed = os.path.join(base_out, "fastqc_trimmed")
        multiqc_trimmed = os.path.join(base_out, "multiqc_trimmed")
        normalized_coverage = os.path.join(base_out, "normalized_coverage")
        peak_files = os.path.join(base_out, "peak_files")
        os.makedirs(fastqc_raw, exist_ok=True)
        os.makedirs(multiqc_raw, exist_ok=True)
        os.makedirs(trimmed_data, exist_ok=True)
        os.makedirs(bam_output, exist_ok=True)
        os.makedirs(fastqc_trimmed, exist_ok=True)
        os.makedirs(multiqc_trimmed, exist_ok=True)
        os.makedirs(normalized_coverage, exist_ok=True)
        os.makedirs(peak_files, exist_ok=True)

        if not self.params["step5"].get("peak_files_dir"):
            self.params["step5"]["peak_files_dir"] = peak_files

        threading.Thread(target=self.execute_pipeline_sequentially, daemon=True).start()

    def execute_pipeline_sequentially(self):
        # Get parameters needed for commands
        base_out = self.params["step1"]["base_output_dir"]
        raw_data = self.params["step1"]["raw_data_dir"]
        samples = self.params["step2"].get("sample_names", [])

        # Define directories
        fastqc_raw = os.path.join(base_out, "fastqc_raw")
        multiqc_raw = os.path.join(base_out, "multiqc_raw")
        trimmed_data = os.path.join(base_out, "trimmed_data")
        bam_output = os.path.join(base_out, "bam_output")
        fastqc_trimmed = os.path.join(base_out, "fastqc_trimmed")
        multiqc_trimmed = os.path.join(base_out, "multiqc_trimmed")
        normalized_coverage = os.path.join(base_out, "normalized_coverage")
        peak_files = os.path.join(base_out, "peak_files")
        threads = self.params.get("threads", 8)

        try:

            # Get selected samples from Step 1
            selected_samples = self.params["step1"].get("selected_samples")

            if selected_samples:
                selected_names = {
                    os.path.basename(f).split("_")[0] for f in selected_samples
                }
                samples = [s for s in self.params["step1"]["auto_sample_names"] if s in selected_names]
            else:
                samples = self.params["step1"]["auto_sample_names"]

            # Step 1: FastQC on raw data
            self.update_output_gui("Checking Step 1: FastQC raw data...\n")
            raw_files = []
            if self.params["step1"].get("selected_samples"):
                raw_files = self.params["step1"]["selected_samples"]
                samples = self.params["step1"]["auto_sample_names"]  # From Step 1
            else:
                # Find all FASTQ files with both extensions
                raw_files = glob.glob(os.path.join(raw_data, "*.fastq.gz")) + glob.glob(
                    os.path.join(raw_data, "*.fq.gz"))

            all_raw_reports_exist = all(
                os.path.exists(os.path.join(fastqc_raw,
                                            os.path.basename(f).replace('.fastq.gz', '_fastqc.html').replace('.fq.gz',
                                                                                                             '_fastqc.html')))
                for f in raw_files
            )

            if not all_raw_reports_exist:
                self.update_output_gui("Running Step 1: FastQC on raw data...\n")
                # Join selected files into a space-separated string
                files_to_process = " ".join(raw_files)
                cmd = f"fastqc -t {threads} {files_to_process} -o {fastqc_raw}"
                if not self.run_blocking_command(cmd):
                    return
            else:
                self.update_output_gui("Skipping Step 1: FastQC raw data (reports exist)\n")

            # Step 2: MultiQC on raw data
            self.update_output_gui("Checking Step 2: MultiQC raw data...\n")
            multiqc_report = os.path.join(multiqc_raw, "multiqc_report.html")
            if not os.path.exists(multiqc_report):
                self.update_output_gui("Running Step 2: MultiQC on raw data...\n")
                cmd = f"multiqc {fastqc_raw} -o {multiqc_raw}"
                if not self.run_blocking_command(cmd):
                    return
            else:
                self.update_output_gui("Skipping Step 2: MultiQC raw data (report exists)\n")

            selected_samples = self.params["step1"].get("selected_samples")
            if not selected_samples:
                raw_data_dir = self.params["step1"]["raw_data_dir"]
                selected_samples = glob.glob(os.path.join(raw_data_dir, "*.fastq.gz")) + \
                                   glob.glob(os.path.join(raw_data_dir, "*.fq.gz"))

            # Step 3: Trimming
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 3: Trim Galore...\n")
                for sample in samples:
                    output_r1 = os.path.join(trimmed_data, f"{sample}_val_1.fq.gz")
                    output_r2 = os.path.join(trimmed_data, f"{sample}_val_2.fq.gz")

                    if not (os.path.exists(output_r1) and os.path.exists(output_r2)):
                        sample_files = [f for f in selected_samples if sample in os.path.basename(f)]

                        r1_list = [f for f in sample_files if
                                   re.search(r"_R?1(?:_001)?\.f(ast)?q\.gz$", os.path.basename(f))]
                        r2_list = [f for f in sample_files if
                                   re.search(r"_R?2(?:_001)?\.f(ast)?q\.gz$", os.path.basename(f))]

                        if not r1_list or not r2_list:
                            self.update_output_gui(f"Warning: Could not find both pairs for sample {sample}\n")
                            continue

                        r1 = r1_list[0]
                        r2 = r2_list[0]

                        self.update_output_gui(f"Running Trim Galore for {sample}...\n")
                        cmd = (f"trim_galore --basename {sample} --gzip  --cores {threads} --paired "
                               f"{r1} {r2} --output_dir {trimmed_data}")
                        if not self.run_blocking_command(cmd):
                            return
                    else:
                        self.update_output_gui(f"Skipping Trim Galore for {sample} (trimmed files exist)\n")
            else:
                self.update_output_gui("Skipping Step 3: Trim Galore (User chose to use raw data)\n")

            # Step 4: FastQC on trimmed data
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 4: FastQC trimmed data...\n")
                trimmed_files = []
                for sample in samples:
                    trimmed_files.append(os.path.join(trimmed_data, f"{sample}_val_1.fq.gz"))
                    trimmed_files.append(os.path.join(trimmed_data, f"{sample}_val_2.fq.gz"))

                all_trimmed_reports_exist = all(
                    os.path.exists(os.path.join(fastqc_trimmed, os.path.basename(f).replace('.fq.gz', '_fastqc.html')))
                    for f in trimmed_files
                )
                if not all_trimmed_reports_exist:
                    self.update_output_gui("Running Step 4: FastQC on trimmed data...\n")
                    cmd = f"fastqc -t {threads} {trimmed_data}/*.fq.gz -o {fastqc_trimmed}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui("Skipping Step 4: FastQC trimmed data (reports exist)\n")
            else:
                self.update_output_gui(
                    "Skipping Steps 4 & 5: FastQC and MultiQC on trimmed data (trimming was skipped).\n")

            # Step 5: MultiQC on trimmed data
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 5: MultiQC trimmed data...\n")
                multiqc_report = os.path.join(multiqc_trimmed, "multiqc_report.html")
                if not os.path.exists(multiqc_report):
                    self.update_output_gui("Running Step 5: MultiQC on trimmed data...\n")
                    cmd = f"multiqc {fastqc_trimmed} -o {multiqc_trimmed}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui("Skipping Step 5: MultiQC trimmed data (report exists)\n")

            # ====== Ensembl Reference Download and Index Build ====== #
            ensembl_release = "111"
            sp = self.params["step3"]["ensembl_species"]
            cap = self.params["step3"]["ensembl_cap"]
            assembly = self.params["step3"]["ensembl_assembly"]
            ref_dir = self.params["step3"]["ref_dir"]
            bt2_base = self.params["step3"]["bt2_base"]
            genome_version = self.params["step3"]["genome_version"]
            fa_file = os.path.join(ref_dir, f"{genome_version}.fa")
            fa_gz = f"{fa_file}.gz"
            index_extensions = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

            if not all(os.path.exists(f"{bt2_base}.{ext}") for ext in index_extensions):
                self.update_output_gui(f"Processing reference genome: {genome_version}\n")

                # Create reference directory if needed
                os.makedirs(ref_dir, exist_ok=True)

                # Check for existing files
                if not os.path.exists(fa_file):
                    if os.path.exists(fa_gz):
                        self.update_output_gui("Found compressed FASTA, decompressing...\n")
                        subprocess.run(f"gzip -d {fa_gz}", shell=True, check=True)
                    else:
                        # Download with aria2c (resume support)
                        self.update_output_gui("Downloading reference genome with aria2c...\n")
                        url = (f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/"
                               f"{sp}/dna/{cap}.{assembly}.dna.toplevel.fa.gz")

                        download_cmd = (
                            f"aria2c -c -x 16 -s 16 -d {ref_dir} "
                            f"-o {genome_version}.fa.gz {url}"
                        )
                        try:
                            subprocess.run(download_cmd, shell=True, check=True)
                            subprocess.run(f"gzip -d {fa_gz}", shell=True, check=True)
                        except subprocess.CalledProcessError as e:
                            self.show_error_gui(f"Download failed: {str(e)}")
                            return

                # Build index if FASTA exists
                if os.path.exists(fa_file):
                    self.update_output_gui("Building Bowtie2 index...\n")
                    try:
                        subprocess.run(
                            f"bowtie2-build --threads {threads} {fa_file} {bt2_base}",
                            shell=True,
                            check=True
                        )
                    except subprocess.CalledProcessError as e:
                        self.show_error_gui(f"Index build failed: {str(e)}")
                        return
                else:
                    self.show_error_gui("FASTA file missing after download/decompression")
                    return
            else:
                self.update_output_gui("Using existing reference index\n")

            # ======== Download GTF ========
            gtf_url = (
                f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}/gtf/"
                f"{sp}/{cap}.{assembly}.{ensembl_release}.gtf.gz"
            )
            gtf_gz = os.path.join(ref_dir, f"{genome_version}.gtf.gz")
            gtf_file = os.path.join(ref_dir, f"{genome_version}.gtf")

            if not os.path.exists(gtf_file):
                if os.path.exists(gtf_gz):
                    self.update_output_gui("Found compressed GTF, decompressing...\n")
                    subprocess.run(f"gzip -d {gtf_gz}", shell=True, check=True)
                else:
                    self.update_output_gui("Downloading GTF...\n")
                    download_cmd = (
                        f"aria2c -c -x 16 -s 16 -d {ref_dir} "
                        f"-o {genome_version}.gtf.gz {gtf_url}"
                    )
                    try:
                        subprocess.run(download_cmd, shell=True, check=True)
                        subprocess.run(f"gzip -d {gtf_gz}", shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        self.show_error_gui(f"GTF download failed: {str(e)}")
                        return
            else:
                self.update_output_gui("GTF file already exists.\n")

            # ====== Obtain TSS file ========= (associated helper function is present after this method)
            tss_bed = os.path.join(ref_dir, f"{genome_version}_TSS.bed")
            if not os.path.exists(tss_bed):
                self.update_output_gui("Generating TSS BED file from GTF...\n")
                try:
                    self.gtf_to_tss_bed(
                        gtf_file=gtf_file,
                        output_bed=tss_bed,
                        genome_assembly=genome_version,
                        gene_biotype="protein_coding"
                    )
                except Exception as e:
                    self.show_error_gui(f"Failed to create TSS BED file: {str(e)}")
                    return

            # Step 6: Alignment and Sorting
            self.update_output_gui("Checking Step 6: Alignment...\n")

            selected_samples = self.params["step1"].get("selected_samples")
            auto_sample_names = self.params["step1"].get("auto_sample_names", [])
            sample_file_map = {}

            if selected_samples:
                # User selected specific FASTQ files
                for r1_path in selected_samples:
                    base_r1 = os.path.basename(r1_path)

                    # Match R1 file and extract sample name
                    match = re.match(r"^(.*?)(_R?1(?:_001)?|_1)\.f(?:ast)?q\.gz$", base_r1)
                    if not match:
                        continue  # Skip unrecognized files

                    sample_name = match.group(1)
                    r2_path = re.sub(r'(_R?)1(\d*)', r'\g<1>2\g<2>', r1_path)

                    if self.params["step2"].get("skip_trimming", False):
                        sample_file_map[sample_name] = {"r1": r1_path, "r2": r2_path}
                    else:
                        trimmed_r1 = os.path.join(trimmed_data, f"{sample_name}_val_1.fq.gz")
                        trimmed_r2 = os.path.join(trimmed_data, f"{sample_name}_val_2.fq.gz")

                        sample_file_map[sample_name] = {"r1": trimmed_r1, "r2": trimmed_r2}

                samples = list(sample_file_map.keys())

            else:
                samples = auto_sample_names
                for sample in samples:
                    if self.params["step2"].get("skip_trimming", False):
                        r1_candidates = glob.glob(os.path.join(raw_data, f"{sample}_*R1*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*1*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*R1*.fq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*1*.fq.gz"))

                        r2_candidates = glob.glob(os.path.join(raw_data, f"{sample}_*R2*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*2*.fastq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*R2*.fq.gz")) + \
                                        glob.glob(os.path.join(raw_data, f"{sample}_*2*.fq.gz"))

                        if not r1_candidates or not r2_candidates:
                            raise ValueError(f"Could not find paired files for sample {sample} in {raw_data}")

                        r1 = sorted(r1_candidates)[0]
                        r2 = sorted(r2_candidates)[0]

                    else:
                        r1 = os.path.join(trimmed_data, f"{sample}_val_1.fq.gz")
                        r2 = os.path.join(trimmed_data, f"{sample}_val_2.fq.gz")

                    sample_file_map[sample] = {"r1": r1, "r2": r2}

            # Main alignment loop
            for sample in samples:
                coord_sorted_bam = os.path.join(bam_output, f"{sample}.sort.bam")
                name_sorted_bam = os.path.join(bam_output, f"{sample}.name.sorted.bam")

                if not (os.path.exists(coord_sorted_bam) and os.path.exists(coord_sorted_bam + ".bai")):
                    self.update_output_gui(f"Aligning {sample}...\n")

                    genome_version = self.params["step3"].get("genome_version", "").strip()

                    r1 = sample_file_map[sample]["r1"]
                    r2 = sample_file_map[sample]["r2"]

                    cmd = (f"bowtie2 -p {threads} --very-sensitive -X 2000 -x {bt2_base} "
                           f"-1 {r1} -2 {r2} "
                           f"| samtools view --threads {threads} -bS - > {bam_output}/{sample}.bam")
                    if not self.run_blocking_command(cmd):
                        return

                    if not os.path.exists(name_sorted_bam):
                        self.update_output_gui(f"Name sorting {sample}...\n")
                        cmd = f"samtools sort -n --threads {threads} -o {name_sorted_bam} {bam_output}/{sample}.bam && rm {bam_output}/{sample}.bam"
                        if not self.run_blocking_command(cmd):
                            return
                    else:
                        self.update_output_gui(f"Skipping name sorting for {sample} (Name-sorted BAM exists)\n")

                    self.update_output_gui(f"Coordinate sorting {sample}...\n")
                    cmd = f"samtools sort --threads {threads} -o {coord_sorted_bam} {name_sorted_bam} && samtools index -@ {threads} {coord_sorted_bam}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui(f"Skipping alignment for {sample} (BAM exists)\n")

            # Paths
            flagstat_path = os.path.join(bam_output, "flagstat_results.txt")
            if os.path.exists(flagstat_path):
                os.remove(flagstat_path)

            bam_files = sorted(glob.glob(os.path.join(bam_output, "*.sort.bam")))

            # Output file header
            with open(flagstat_path, "w") as stats_file:
                stats_file.write(
                    "Sample\tTotal Reads\tMapped (%)\tProperly Paired (%)\tSingletons (%)\tDiscordant (mapQ>=5)\n")

                for bam in bam_files:
                    name = os.path.basename(bam)
                    self.update_output_gui(f"Generating alignment stats for {name}...\n")

                    # Run samtools and capture output
                    result = subprocess.run(
                        ["samtools", "flagstat", bam],
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    output = result.stdout

                    # Extract numbers using regex
                    total_reads = re.search(r"(\d+) \+ \d+ in total", output)
                    mapped = re.search(r"(\d+) \+ \d+ mapped \(([\d\.]+)%", output)
                    properly_paired = re.search(r"(\d+) \+ \d+ properly paired \(([\d\.]+)%", output)
                    singletons = re.search(r"(\d+) \+ \d+ singletons \(([\d\.]+)%", output)
                    discordant = re.search(r"(\d+) \+ \d+ with mate mapped to a different chr \(mapQ>=5\)", output)

                    # Extract and fallback
                    total = total_reads.group(1) if total_reads else "NA"
                    mapped_pct = mapped.group(2) if mapped else "NA"
                    proper_pct = properly_paired.group(2) if properly_paired else "NA"
                    single_pct = singletons.group(2) if singletons else "NA"
                    discordant_reads = discordant.group(1) if discordant else "NA"

                    # Write summary line
                    stats_file.write(f"{name}\t{total}\t{mapped_pct}\t{proper_pct}\t{single_pct}\t{discordant_reads}\n")

            # Step 7: Coverage
            self.update_output_gui("Checking Step 7: Coverage...\n")
            for sample in samples:
                coverage_bw = os.path.join(normalized_coverage, f"{sample}.normalized.bw")
                if not os.path.exists(coverage_bw):
                    self.update_output_gui(f"Generating coverage for {sample}...\n")
                    cmd = f"bamCoverage -p {threads} -of bigwig --normalizeUsing=RPKM -b {bam_output}/{sample}.sort.bam -o {coverage_bw}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui(f"Skipping coverage for {sample} (file exists)\n")

            # Create output directory for coverage profiles (associated helper function after this method)
            coverage_profile_dir = os.path.join(normalized_coverage, "coverage_profiles")
            os.makedirs(coverage_profile_dir, exist_ok=True)

            # expected output file paths
            tss_prefix = os.path.join(coverage_profile_dir, "TSS")
            matrix_file = f"{tss_prefix}_matrix.gz"
            heatmap_file = f"{tss_prefix}_heatmap.pdf"

            # Check if both output files already exist
            if os.path.exists(matrix_file) and os.path.exists(heatmap_file):
                self.update_output_gui("TSS coverage profile and heatmap already exist. Skipping.\n")
            else:
                self.update_output_gui("Generating TSS coverage profiles...\n")
                try:
                    self.run_tss_matrix_and_heatmap(
                        tss_bed=tss_bed,
                        bw_dir=normalized_coverage,
                        out_prefix=tss_prefix
                    )
                except Exception as e:
                    self.show_error_gui(f"Failed to generate TSS matrix and heatmap: {str(e)}")
                    return

            # =================== Step 8: Peak Calling (Genrich or MACS3) =================== #

            # Genrich
            if self.params["step4"].get("peak_caller") == "Genrich":
                bam_dir = bam_output
                treated = self.params["step4"].get("treated_samples", [])
                control = self.params["step4"].get("control_samples", [])
                blacklist = ",".join(self.params["step4"].get("blacklist_files", []))
                # For exclude_chr, if non-empty then join using comma; otherwise, omit the parameter
                exclude_chr_list = self.params["step4"].get("exclude_chr", [])
                if exclude_chr_list:
                    # Create Genrich argument: -e expects comma-separated list.
                    exclude_chr_arg = f"-e {','.join(exclude_chr_list)}"
                else:
                    exclude_chr_arg = ""
                expand_cut = self.params["step4"].get("expand_cut_sites", 100)
                max_q = self.params["step4"].get("max_q_value", 0.05)

                # For blacklist files, include -E flag only when blacklist is provided.
                if blacklist:
                    blacklist_arg = f"-E {blacklist}"
                else:
                    blacklist_arg = ""

                treated_files = ",".join([f"{bam_dir}/{s}.name.sorted.bam" for s in treated])
                control_flag = f"-c {','.join([f'{bam_dir}/{s}.name.sorted.bam' for s in control])}" if control else ""

                if self.params["step4"].get("peak_type") == "merged":
                    output_peak = f"{peak_files}/{self.params['step4'].get('merged_output_name')}.genrich.peak"
                    output_bg = f"{peak_files}/{self.params['step4'].get('merged_output_name')}.genrich.peak.bg"
                    cmd_genrich = (f"Genrich -t {treated_files} {control_flag} -o {output_peak} -f {output_bg} "
                                   f"-v -j -d {expand_cut} -r -q {max_q} {exclude_chr_arg} {blacklist_arg}")
                    self.update_output_gui("Running Step 8: Genrich Peak Calling...\n")
                    if not self.run_blocking_command(cmd_genrich):
                        return
                    self.run_peak_annotation()

                else:
                    for sample in treated:
                        output_peak = f"{peak_files}/{sample}.genrich.peak"
                        output_bg = f"{peak_files}/{sample}.genrich.peak.bg"
                        cmd_genrich = (f"Genrich -t {bam_dir}/{sample}.name.sorted.bam {control_flag} "
                                       f"-o {output_peak} -f {output_bg} -v -j -d {expand_cut} -r "
                                       f"-q {max_q} {exclude_chr_arg} {blacklist_arg}")
                        self.update_output_gui(f"Running Genrich for {sample}...\n")
                        if not self.run_blocking_command(cmd_genrich):
                            return
                    self.run_peak_annotation()

            # MACS3
            else:
                bam_dir = bam_output
                treated = self.params["step4"].get("treated_samples", [])
                control = self.params["step4"].get("control_samples", [])
                blacklist_files = self.params["step4"].get("blacklist_files", [])
                exclude_chr = self.params["step4"].get("exclude_chr", [])
                max_q = self.params["step4"].get("max_q_value", 0.05)
                genome_size = self.params["step4"].get("genome", "hs")

                files_to_delete = set()
                control_filtered_map = {}
                macs3_commands = []

                # Helper to generate BAM processing command
                def process_bam(sample, files_to_delete):
                    name_sorted_bam = os.path.join(bam_dir, f"{sample}.name.sorted.bam")
                    fixed_bam = os.path.join(bam_dir, f"{sample}.fixmate.bam")
                    coord_bam = os.path.join(bam_dir, f"{sample}.coord.bam")
                    dedup_bam = os.path.join(bam_dir, f"{sample}.dedup.bam")
                    filtered_bam = os.path.join(bam_dir, f"{sample}.filtered.dedup.bam")

                    cmd = (
                        f"samtools fixmate -m {name_sorted_bam} {fixed_bam} && "
                        f"samtools sort --threads {threads} -o {coord_bam} {fixed_bam} && "
                        f"samtools markdup -r --threads {threads} {coord_bam} {dedup_bam} && "
                        f"samtools index -@ {threads} {dedup_bam}"
                    )
                    files_to_delete.update([fixed_bam, coord_bam, dedup_bam])
                    return cmd, dedup_bam, filtered_bam

                # Process controls
                for ctrl in control:
                    self.update_output_gui(f"Processing control sample {ctrl}...\n")
                    cmd, dedup, filtered = process_bam(ctrl, files_to_delete)
                    if not self.run_blocking_command(cmd): return
                    try:
                        self.filter_bam(dedup, filtered, blacklist_files, exclude_chr)
                        files_to_delete.add(filtered)
                        control_filtered_map[ctrl] = filtered
                    except Exception as e:
                        self.show_error_gui(f"BAM filtering failed for control sample {ctrl}: {e}")
                        return

                # Process treated + construct MACS3 commands
                for s in treated:
                    self.update_output_gui(f"Processing test sample {s}...\n")
                    cmd, dedup, filtered = process_bam(s, files_to_delete)
                    if not self.run_blocking_command(cmd): return
                    try:
                        self.filter_bam(dedup, filtered, blacklist_files, exclude_chr)
                        files_to_delete.add(filtered)
                    except Exception as e:
                        self.show_error_gui(f"BAM filtering failed for sample {s}: {e}")
                        return

                    if control:
                        for ctrl in control:
                            control_bam = control_filtered_map[ctrl]
                            prefix = os.path.join(peak_files, f"{s}_{ctrl}.macs3.peak")
                            macs_cmd = (
                                f"macs3 callpeak -f BAMPE --call-summits -t {filtered} -c {control_bam} "
                                f"-g {genome_size} -n {prefix} -B -q {max_q}"
                            )
                            macs3_commands.append((f"{s} vs {ctrl}", macs_cmd))
                    else:
                        prefix = os.path.join(peak_files, f"{s}.macs3.peak")
                        macs_cmd = (
                            f"macs3 callpeak -f BAMPE --call-summits -t {filtered} "
                            f"-g {genome_size} -n {prefix} -B -q {max_q}"
                        )
                        macs3_commands.append((s, macs_cmd))

                # Run MACS3
                self.update_output_gui("Running Step 8: MACS3 Peak Calling...\n")
                self.run_macs3_commands(
                    macs3_commands,
                    on_complete=lambda: self.cleanup_and_continue(files_to_delete)
                )

        except Exception as e:
            self.show_error_gui(f"Pipeline failed: {str(e)}")

    def run_macs3_commands(self, commands, on_complete=None):
        def worker():
            for label, cmd in commands:
                self.update_output_gui(f"\n Calling peaks for: {label}\n{cmd}\n")
                try:
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                            text=True)
                    for line in proc.stdout:
                        self.root.after(0, lambda l=line: self.update_output_gui(l))
                    exit_code = proc.wait()
                    if exit_code != 0:
                        self.root.after(0, lambda: self.show_error_gui(f"MACS3 failed for {label}"))
                        return
                except Exception as e:
                    self.root.after(0, lambda: self.show_error_gui(f"MACS3 crashed for {label}: {e}"))
                    return
            self.root.after(0, lambda: self.update_output_gui("✅ MACS3 Peak Calling completed.\n"))
            if on_complete:
                self.root.after(0, on_complete)

        threading.Thread(target=worker, daemon=True).start()

    def cleanup_and_continue(self, files_to_delete):
        self.update_output_gui("Cleaning up temporary BAM files...\n")
        for f in files_to_delete:
            try:
                os.remove(f)
            except Exception as e:
                self.update_output_gui(f"Warning: Could not delete file {f}: {e}\n")
        self.run_peak_annotation()

    def run_peak_annotation(self):
        r_script_path = resource_filename('chromacs', 'annotate_peaks_5.R')
        if not os.path.exists(r_script_path):
            self.show_error_gui(f"R script 'annotate_peaks_5.R' not found at:\n{r_script_path}")
            return

        genome_version = self.params["step3"]["genome_version"]
        ref_dir = self.params["step3"]["ref_dir"]
        peak_files = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        cmd = f"Rscript {r_script_path} {peak_files} {genome_version} {ref_dir}"
        self.update_output_gui(f"Running Peak Annotation with {r_script_path}...\n")

        def worker():
            try:
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True)
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                proc.wait()
                self.root.after(0, lambda: self.diffbind_btn.config(state=tk.NORMAL))
                timestamp = self.get_timestamp()
                self.update_output_gui(f"\n✅ Pipeline execution complete! [{timestamp}]\n")
                messagebox.showinfo("Pipeline Finished", "✅ All steps have completed successfully!")
            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(f"Peak annotation failed: {e}"))

        threading.Thread(target=worker, daemon=True).start()

    # this is required for MACS3 filtering
    def filter_bam(self, input_bam, output_bam, blacklist_files, exclude_chr):
        import tempfile
        import shutil

        tmp_dir = tempfile.mkdtemp()
        step1 = os.path.join(tmp_dir, "step1.bam")

        try:
            # Step 1: Blacklist filtering
            if blacklist_files:
                blacklist_str = " ".join(blacklist_files)
                cmd1 = f"bedtools intersect -v -abam {input_bam} -b {blacklist_str} > {step1}"
            else:
                shutil.copy(input_bam, step1)
                cmd1 = None

            if cmd1:
                self.update_output_gui(f"Running blacklist filtering:\n{cmd1}\n")
                subprocess.run(cmd1, shell=True, check=True)

            # Step 2: Exclude chromosomes
            if exclude_chr:
                pat = "|".join(re.escape(c) for c in exclude_chr)
                cmd2 = (
                    f"samtools view -h {step1} | "
                    f"awk '$1 ~ /^@/ || $3 !~ /^({pat})$/' | "
                    f"samtools view -b - > {output_bam}"
                )
            else:
                shutil.copy(step1, output_bam)
                cmd2 = None

            if cmd2:
                self.update_output_gui(f"Running exclude chromosome filter:\n{cmd2}\n")
                subprocess.run(cmd2, shell=True, check=True)

        except subprocess.CalledProcessError as e:
            raise Exception(f"BAM filtering failed: {str(e)}")
        finally:
            shutil.rmtree(tmp_dir)

    # required for coverage matrix calculation and plotting of bigwigs
    def gtf_to_tss_bed(
            self,
            gtf_file: str,
            output_bed: str,
            genome_assembly: str,
            gene_biotype: str = "protein_coding"
    ):
        """
        Converts GTF to TSS BED, excluding mitochondrial chromosomes dynamically based on genome assembly.

        Args:
            gtf_file: Input GTF path (can be gzipped).
            output_bed: Output BED path.
            genome_assembly: Genome assembly name like "GRCh38", "GRCm39".
            gene_biotype: Filter gene biotype.
        """

        mito_chromosomes = {
            "GRCh38": {"chrM", "MT", "M"},
            "GRCm39": {"chrM", "MT", "M"},
            "mRatBN7.2": {"chrM", "MT", "M"},
            "ARS-UCD1.3": {"chrM", "MT", "M"},  # Bos taurus mitochondrion is usually "MT" or "chrM"
            "Sscrofa11.1": {"chrM", "MT", "M"},
            "GRCg7b": {"chrM", "MT", "M"},  # Gallus gallus (chicken) mito often "chrM" or "MT"
            "Pan_tro_3.0": {"chrM", "MT", "M"},  # Chimpanzee mitochondrion likely similar to human
            "ROS_Cfam_1.0": {"chrM", "MT", "M"},  # Dog mitochondrion
            "ARS1": {"chrM", "MT", "M"},  # Capra hircus (goat) mitochondrion, check assembly docs if needed
            "CVASU_BBG_1.0": {"chrM", "MT", "M"},  # Another goat assembly, same as above
            "OryCun2.0": {"chrM", "MT", "M"},  # Rabbit mitochondrion
            "gorGor4": {"chrM", "MT", "M"},  # Gorilla mitochondrion
            "Mmul_10": {"chrM", "MT", "M"},  # Macaca mulatta (rhesus macaque)

            # Fish & Amphibians
            "GRCz11": {"chrM", "MT", "M"},
            "UCB_Xtro_10.0": {"chrM", "MT", "M"},  # Xenopus tropicalis
            "Ssal_v3.1": {"chrM", "MT", "M"},  # Salmo salar (Atlantic salmon)

            # Invertebrates
            "BDGP6": {"chrM", "mitochondrion", "MT"},
            "WBcel235": {"MtDNA", "chrM"},
        }

        # Keyword-based blacklist (safe for most organisms)
        non_primary_keywords = {"random", "alt", "hap", "fix", "patch", "scaffold", "gl", "un", "ki"}

        def is_primary_contig(chrom: str) -> bool:
            chrom = chrom.lower()
            return not any(k in chrom for k in non_primary_keywords)

        exclude_chroms = mito_chromosomes.get(genome_assembly, set())

        with (
                gzip.open(gtf_file, "rt") if gtf_file.endswith(".gz") else open(gtf_file, "r")
        ) as infile, open(output_bed, "w") as outfile:

            seen_genes = set()

            for line in infile:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9 or cols[2].lower() != "gene":
                    continue

                chrom = cols[0]

                # Exclude mitochondrial chromosomes for this genome assembly
                if chrom in exclude_chroms:
                    continue

                # Optionally skip non-primary contigs (unplaced scaffolds, alt loci)
                if not is_primary_contig(chrom):
                    continue

                # Parse attributes
                attrs = {}
                for pair in cols[8].strip().split(";"):
                    pair = pair.strip()
                    if not pair:
                        continue
                    if "\"" in pair:
                        key, val = pair.split(" ", 1)
                        val = val.strip("\"")
                    else:
                        key, val = pair.split("=") if "=" in pair else (pair, "")
                    attrs[key.strip()] = val.strip()

                # Filter by gene biotype
                biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or attrs.get("biotype")
                if not biotype or biotype.lower() != gene_biotype.lower():
                    continue

                strand = cols[6]
                start = int(cols[3])
                end = int(cols[4])
                tss = start if strand == "+" else end

                bed_start = tss - 1
                bed_end = bed_start + 1

                gene_id = attrs.get("gene_id") or attrs.get("ID")
                if not gene_id or gene_id in seen_genes:
                    continue
                seen_genes.add(gene_id)

                outfile.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\t.\t{strand}\n")

    def run_tss_matrix_and_heatmap(self, tss_bed: str, bw_dir: str, out_prefix: str):

        bw_files = glob.glob(os.path.join(bw_dir, "*.normalized.bw"))
        if not bw_files:
            raise RuntimeError(f"No BigWig files found in {bw_dir}")

        matrix_file = f"{out_prefix}_matrix.gz"
        cmd_mat = [
                      "computeMatrix", "reference-point",
                      "--referencePoint", "center",
                      "-b", "3000", "-a", "3000",
                      "-R", tss_bed,
                      "-S"
                  ] + bw_files + [
                      "--skipZeros",
                      "-p", "max",
                      "-o", matrix_file
                  ]
        subprocess.run(cmd_mat, check=True)

        heatmap_file = f"{out_prefix}_heatmap.pdf"
        cmd_heat = [
            "plotHeatmap",
            "-m", matrix_file,
            "-out", heatmap_file,
            "--colorMap", "Reds",
            "--heatmapHeight", "15",
            "--heatmapWidth", "20",
            "--xAxisLabel", "Distance from TSS (bp)",
            "--refPointLabel", "TSS",
            "--averageTypeSummaryPlot", "mean",
            "--plotType", "se",
            "--legendLocation", "best"
        ]
        subprocess.run(cmd_heat, check=True)

        return matrix_file, heatmap_file

    # =================== DiffBind Analysis Logic =================== #

    def launch_diffbind_config(self):
        self.diffbind_window = tk.Toplevel(self.root)
        self.diffbind_window.title("DiffBind Configuration")
        self.diffbind_window.geometry("1000x700")

        # Get all samples for ControlID dropdown
        self.all_samples = self.params["step2"].get("sample_names", [])

        # Peak file selection
        tk.Label(self.diffbind_window, text="1. Select Peak Files:",
                 font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, padx=10, pady=5)
        self.peak_listbox = tk.Listbox(self.diffbind_window, selectmode=tk.MULTIPLE, exportselection=False, height=10,
                                       width=80)
        self.peak_listbox.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.diffbind_window, text="Browse Peaks",
                  command=self.populate_peak_list, bg='gray').grid(row=0, column=2)

        # Configuration grid
        self.config_frame = tk.Frame(self.diffbind_window)
        self.config_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=10, sticky='nw')

        # Instruction label (centered)
        instruction_1 = (
            "Note: You can skip this Control section, if there is none\n"
            "Please keep in mind that these are not biological controls, rather technical controls\n"
            "Do NOT put untreated condition as control here\n"
        )

        tk.Label(
            self.diffbind_window, text=instruction_1,
            wraplength=650, justify=tk.LEFT, anchor="w",
            font=(self.roboto_font, 9)
        ).grid(row=1, column=2, padx=10, sticky="w")

        # FDR Threshold Input
        tk.Label(self.diffbind_window, text="FDR Threshold (0–1):", font=(self.roboto_font, 10)).grid(row=2, column=0,
                                                                                                      padx=10,
                                                                                                      sticky='w')
        self.fdr_threshold = tk.DoubleVar(value=0.05)
        tk.Spinbox(self.diffbind_window, from_=0.0, to=1.0, increment=0.01, textvariable=self.fdr_threshold,
                   width=5).grid(row=2, column=1, sticky='w')

        # Thread count input
        tk.Label(self.diffbind_window, text="Number of Threads (default = 8):", font=(self.roboto_font, 10)).grid(
            row=4, column=0, sticky='w', padx=10)
        self.diffbind_threads = tk.StringVar(value="8")  # default is 8
        tk.Entry(self.diffbind_window, textvariable=self.diffbind_threads, width=5).grid(
            row=4, column=1, sticky='w')

        # Instruction label (centered)
        instruction_1 = (
            "Note: DiffBind is very memory intensive tool, so adjust your thread here accordingly\n"
            "      If run fails with error in the core setup, adjust (decrease thread count) and rerun DiffBind\n"
        )
        tk.Label(
            self.diffbind_window, text=instruction_1,
            wraplength=650, justify=tk.LEFT, anchor="w",
            font=(self.roboto_font, 9)
        ).grid(row=6, column=1, padx=10, sticky="w")

        # Run button
        tk.Button(self.diffbind_window, text="Run DiffBind",
                  command=self.validate_and_run_diffbind, bg="yellow green").grid(row=8, column=1, pady=10)

        # Initialize metadata stores
        self.conditions = {}
        self.replicates = {}
        self.controls = {}

        self.peak_listbox.bind('<<ListboxSelect>>', lambda e: self.build_param_rows())

        # annotate diffbind results
        tk.Button(
            self.diffbind_window,
            text="Annotate Results",
            command=self.launch_diffbind_annotation,
            bg="cyan"
        ).grid(row=9, column=1, pady=10)

    def populate_peak_list(self):
        self.peak_listbox.delete(0, tk.END)
        peak_dir = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        peak_files = glob.glob(os.path.join(peak_dir, "*.narrowPeak")) \
                     + glob.glob(os.path.join(peak_dir, "*.genrich.peak"))
        for pf in sorted(peak_files):
            name = os.path.basename(pf)
            self.peak_listbox.insert(tk.END, name)
            # Default vars
            self.conditions[name] = tk.StringVar(value="treated")
            self.replicates[name] = tk.IntVar(value=1)
            self.controls[name] = tk.StringVar(value="")

    def build_param_rows(self):
        # Clear existing rows
        for child in self.config_frame.winfo_children():
            child.destroy()

        # Header
        headers = ['Peak File', 'Condition', 'Replicate', 'Control']
        for col, h in enumerate(headers):
            tk.Label(self.config_frame, text=h, font=(self.roboto_font, 9, 'bold')).grid(row=0, column=col, padx=5)

        # Rows for each selected peak
        for i, idx in enumerate(self.peak_listbox.curselection(), start=1):
            peak = self.peak_listbox.get(idx)
            # Initialize variables if not exists
            if peak not in self.conditions:
                self.conditions[peak] = tk.StringVar(value="treated")
            if peak not in self.replicates:
                self.replicates[peak] = tk.IntVar(value=1)
            if peak not in self.controls:
                self.controls[peak] = tk.StringVar(value="")

            # Peak file name
            tk.Label(self.config_frame, text=peak).grid(row=i, column=0, sticky='w')

            # Condition dropdown
            ttk.Combobox(self.config_frame, textvariable=self.conditions[peak],
                         values=["treated", "untreated"], state="readonly").grid(row=i, column=1)

            # Replicate spinbox
            tk.Spinbox(self.config_frame, from_=1, to=10,
                       textvariable=self.replicates[peak], width=5).grid(row=i, column=2)

            # Control ID dropdown (all samples)
            ttk.Combobox(self.config_frame, textvariable=self.controls[peak],
                         values=self.all_samples, state="readonly").grid(row=i, column=3)

    def validate_and_run_diffbind(self):
        if not self.peak_listbox.curselection():
            return self.show_error_gui("Please select at least one peak file.")

        has_any_control = False
        missing_controls = []

        # Validate required fields
        for idx in self.peak_listbox.curselection():
            peak = self.peak_listbox.get(idx)
            if not self.conditions[peak].get():
                return self.show_error_gui(f"Condition missing for {peak}")
            if not self.replicates[peak].get():
                return self.show_error_gui(f"Replicate missing for {peak}")
            if self.controls[peak].get():
                has_any_control = True
            else:
                missing_controls.append(peak)

        # Warn the user if some peaks have control and others do not
        if has_any_control and missing_controls:
            peak_list = "\n".join(missing_controls)
            proceed = messagebox.askyesno(
                "Mixed Control Usage",
                f"Some peak files are missing control inputs:\n\n{peak_list}\n\n"
                "Do you want to continue without controls for those?"
            )
            if not proceed:
                return

        fdr = self.fdr_threshold.get()
        if not (0.0 <= fdr <= 1.0):
            return self.show_error_gui("FDR threshold must be between 0 and 1.")

        self.run_diffbind_analysis()

    def generate_metadata(self):
        rows = []
        base = self.params["step1"]["base_output_dir"]
        has_control = False  # Flag to check if any control is used

        for idx in self.peak_listbox.curselection():
            peak_filename = self.peak_listbox.get(idx)
            sample_id_with_suffix = os.path.basename(peak_filename).split('.')[0]

            if 'macs3' in peak_filename.lower():
                sample_id = sample_id_with_suffix.split('_')[0]
                peak_caller = 'narrow'
            else:
                sample_id = sample_id_with_suffix
                peak_caller = 'bed'

            condition = self.conditions[peak_filename].get()
            replicate = self.replicates[peak_filename].get()
            control_id = self.controls[peak_filename].get()
            if control_id: has_control = True

            bam_reads = os.path.join(base, "bam_output", f"{sample_id}.sort.bam")
            peak_path = os.path.join(base, "peak_files", peak_filename)
            bam_control = os.path.join(base, "bam_output", f"{control_id}.sort.bam") if control_id else ""

            row = {
                "SampleID": sample_id,
                "bamReads": bam_reads,
                "Condition": condition,
                "Replicate": replicate,
                "Peaks": peak_path,
                "PeakCaller": peak_caller
            }

            if control_id:
                row["ControlID"] = control_id
                row["bamControl"] = bam_control

            rows.append(row)

        df = pd.DataFrame(rows)

        # Reorder columns if controls are included
        if has_control:
            cols = ["SampleID", "bamReads", "Condition", "Replicate",
                    "Peaks", "PeakCaller", "ControlID", "bamControl"]
        else:
            cols = ["SampleID", "bamReads", "Condition", "Replicate",
                    "Peaks", "PeakCaller"]

        return df[cols]

    def run_diffbind_analysis(self):
        # 1. Prepare output directory
        base_dir = self.params["step1"]["base_output_dir"]
        out_dir = os.path.join(base_dir, "diffbind_results")
        os.makedirs(out_dir, exist_ok=True)

        # 2. Generate and save the metadata sheet
        df = self.generate_metadata()
        out_csv = os.path.join(out_dir, "diffbind_samplesheet.csv")
        df.to_csv(out_csv, index=False)

        # 3. Build the Rscript command
        r_script_path = resource_filename('chromacs', 'diffbind_3.R')
        fdr = str(self.fdr_threshold.get())
        threads = self.diffbind_threads.get().strip()
        threads = threads if threads.isdigit() and int(threads) > 0 else "8"  # Fallback to 8 if invalid

        cmd = ["Rscript", r_script_path, out_csv, out_dir, fdr, threads]

        # 4. Define the worker that actually calls R
        def worker():
            try:
                # Launch Rscript, merge stderr into stdout
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1  # line-buffered
                )

                # Read each line as it comes and post back to the GUI
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))

                exit_code = proc.wait()

                # On success
                if exit_code == 0:
                    def on_success():
                        timestamp = self.get_timestamp()
                        self.update_output_gui(f"✅ DiffBind analysis completed successfully! [{timestamp}]\n")
                        messagebox.showinfo(
                            "DiffBind Completed",
                            f"Results are in:\n{out_dir}"
                        )

                    self.root.after(0, on_success)

                # On failure
                else:
                    def on_error():
                        msg = f" DiffBind failed with exit code {exit_code}\n"
                        self.update_output_gui(msg)
                        self.show_error_gui(msg)

                    self.root.after(0, on_error)

            except Exception as e:
                # Catch unexpected errors
                def on_exception():
                    msg = f" Unexpected error running DiffBind:\n{e}\n"
                    self.update_output_gui(msg)
                    self.show_error_gui(msg)

                self.root.after(0, on_exception)

        # 5. Kick off the worker thread
        threading.Thread(target=worker, daemon=True).start()

    def launch_diffbind_annotation(self):
        base_dir = self.params["step1"]["base_output_dir"]
        diffbind_dir = os.path.join(base_dir, "diffbind_results")
        ref_dir = self.params["step3"]["ref_dir"]
        assembly = self.params["step3"]["genome_version"]

        r_script = resource_filename('chromacs', 'annotate_diffbind.R')

        cmd = [
            "Rscript", r_script,
            diffbind_dir,
            assembly,
            ref_dir
        ]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                proc.wait()
                self.root.after(0, lambda: messagebox.showinfo(
                    "Success", "DiffBind annotation complete!"
                ))
            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(str(e)))

        threading.Thread(target=worker, daemon=True).start()

    # =================== NoiSeq Analysis Logic =================== #

    def launch_noisq_config(self):
        self.noisq_window = tk.Toplevel(self.root)
        self.noisq_window.title("NOISeq Configuration")
        self.noisq_window.geometry("1000x600")

        # Peak file selection
        tk.Label(self.noisq_window, text="Select Peak Files (1 per sample):",
                 font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, padx=10, pady=5)

        self.noisq_peak_listbox = tk.Listbox(self.noisq_window, selectmode=tk.MULTIPLE, exportselection=False,
                                             height=10, width=80)
        self.noisq_peak_listbox.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.noisq_window, text="Browse Peaks",
                  command=self.populate_noisq_peak_list, bg='gray').grid(row=0, column=2)

        # Metadata assignment
        self.noisq_meta_frame = tk.Frame(self.noisq_window)
        self.noisq_meta_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=10)

        self.noisq_conditions = {}

        self.noisq_peak_listbox.bind('<<ListboxSelect>>', lambda e: self.build_noisq_param_rows())

        # NOISeq q-value threshold
        tk.Label(self.noisq_window, text="Confidence Threshold (NOISeq q, 0–1):",
                 font=(self.roboto_font, 10)).grid(row=2, column=0, padx=10, sticky='w')
        self.noisq_qvalue = tk.DoubleVar(value=0.9)
        tk.Spinbox(self.noisq_window, from_=0.0, to=1.0, increment=0.01,
                   textvariable=self.noisq_qvalue, width=5).grid(row=2, column=1, sticky='w')

        # Run button
        tk.Button(self.noisq_window, text="Run NOISeq",
                  command=self.run_noisq_gui, bg="light green").grid(row=3, column=1, pady=10)

        # annotate button
        tk.Button(
            self.noisq_window,
            text="Annotate NOISeq Results",
            command=self.launch_noisq_annotation,
            bg="cyan"
        ).grid(row=4, column=1, pady=10)

    def populate_noisq_peak_list(self):
        self.noisq_peak_listbox.delete(0, tk.END)
        peak_dir = os.path.join(self.params["step1"]["base_output_dir"], "peak_files")
        peak_files = glob.glob(os.path.join(peak_dir, "*.narrowPeak")) \
                     + glob.glob(os.path.join(peak_dir, "*.genrich.peak"))
        for pf in sorted(peak_files):
            name = os.path.basename(pf)
            self.noisq_peak_listbox.insert(tk.END, name)
            self.noisq_conditions[name] = tk.StringVar(value="treated")

    def build_noisq_param_rows(self):
        for child in self.noisq_meta_frame.winfo_children():
            child.destroy()

        tk.Label(self.noisq_meta_frame, text="Peak File").grid(row=0, column=0, padx=5)
        tk.Label(self.noisq_meta_frame, text="Condition").grid(row=0, column=1, padx=5)

        for i, idx in enumerate(self.noisq_peak_listbox.curselection(), start=1):
            peak = self.noisq_peak_listbox.get(idx)
            tk.Label(self.noisq_meta_frame, text=peak).grid(row=i, column=0, sticky='w')
            ttk.Combobox(self.noisq_meta_frame, textvariable=self.noisq_conditions[peak],
                         values=["treated", "untreated"], state="readonly").grid(row=i, column=1)

    def generate_noisq_metadata(self):
        rows = []
        base = self.params["step1"]["base_output_dir"]

        for idx in self.noisq_peak_listbox.curselection():
            peak_filename = self.noisq_peak_listbox.get(idx)
            sample_id_with_suffix = os.path.basename(peak_filename).split('.')[0]

            # Apply same fix as DiffBind
            if 'macs3' in peak_filename.lower():
                # Take only the first SRR/identifier (before underscore)
                sample_id = sample_id_with_suffix.split('_')[0]
            else:
                sample_id = sample_id_with_suffix

            peak_path = os.path.join(base, "peak_files", peak_filename)
            bam_reads = os.path.join(base, "bam_output", f"{sample_id}.sort.bam")

            # You may already have fields like condition/replicate set via UI
            condition = self.noisq_conditions[peak_filename].get()

            rows.append({
                "SampleID": sample_id,
                "Condition": condition,
                "Peaks": peak_path,
                "bamReads": bam_reads,
            })

        return pd.DataFrame(rows)

    def run_noisq_gui(self):
        selected = self.noisq_peak_listbox.curselection()
        if not selected:
            self.show_error_gui("Please select at least one peak file for NOISeq.")
            return

        base_dir = self.params["step1"]["base_output_dir"]
        peak_dir = os.path.join(base_dir, "peak_files")
        out_dir = os.path.join(base_dir, "noisq_results")
        os.makedirs(out_dir, exist_ok=True)

        metadata_df = self.generate_noisq_metadata()
        meta_csv = os.path.join(out_dir, "metadata.csv")
        metadata_df.to_csv(meta_csv, index=False)
        threads = self.params.get("threads", 8)

        # Generate consensus peaks and count matrix
        try:
            # Generate consensus peaks
            all_peaks = [os.path.join(peak_dir, self.noisq_peak_listbox.get(idx)) for idx in selected]
            merged_bed = os.path.join(out_dir, "consensus_peaks.bed")

            # Use absolute paths and check command success
            cat_cmd = f"cat {' '.join(all_peaks)} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > {merged_bed}"
            result = subprocess.run(cat_cmd, shell=True, check=True, capture_output=True, text=True)
            self.update_output_gui(result.stdout)

            # Generate SAF file
            saf_file = os.path.join(out_dir, "consensus_peaks.saf")
            with open(merged_bed, 'r') as infile, open(saf_file, 'w') as outfile:
                outfile.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                for i, line in enumerate(infile, 1):
                    parts = line.strip().split('\t')
                    outfile.write(f"peak{i}\t{parts[0]}\t{parts[1]}\t{parts[2]}\t.\n")

            # Run featureCounts
            bam_files = metadata_df["bamReads"].tolist()
            count_matrix = os.path.join(out_dir, "peak_counts.txt")
            feature_cmd = (
                f"featureCounts -a {saf_file} -F SAF -T {threads} -p -B -C "
                f"-o {count_matrix} {' '.join(bam_files)}"
            )
            result = subprocess.run(feature_cmd, shell=True, check=True, capture_output=True, text=True)
            self.update_output_gui(result.stdout)

            # Process counts matrix
            counts_df = pd.read_csv(count_matrix, sep='\t', comment='#', skiprows=1)
            counts_df = counts_df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
            counts_df = counts_df.rename(columns={"Geneid": "peak_id"})

            # Fix column headers to use SampleIDs instead of full BAM paths
            bam_to_sample = dict(zip(metadata_df["bamReads"], metadata_df["SampleID"]))
            # Leave 'peak_id' column intact, rename only BAM paths
            new_columns = []
            for col in counts_df.columns:
                if col == "peak_id":
                    new_columns.append(col)
                else:
                    new_columns.append(bam_to_sample.get(col, col))
            counts_df.columns = new_columns

            # Save the cleaned count matrix
            peak_matrix_clean = os.path.join(out_dir, "peak_matrix.csv")
            counts_df.to_csv(peak_matrix_clean, index=False)


        except subprocess.CalledProcessError as e:
            self.show_error_gui(f"Command failed: {e.cmd}\nError: {e.stderr}")
            return

        # Run R script
        r_script = resource_filename('chromacs', 'noisq_atac.R')
        output_xlsx = os.path.join(out_dir, "noisq_results.xlsx")
        counts_csv = os.path.join(out_dir, "peak_matrix.csv")
        qval = str(self.noisq_qvalue.get())
        cmd = ["Rscript", r_script, counts_csv, meta_csv, output_xlsx, qval]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                for line in iter(proc.stdout.readline, ''):
                    self.update_output_gui(line)
                proc.wait()
                if proc.returncode == 0:
                    timestamp = self.get_timestamp()
                    self.update_output_gui(f"✅ NOISeq analysis completed! [{timestamp}]\n")
                    messagebox.showinfo("Success", f"Results saved to:\n{output_xlsx}")
                else:
                    self.show_error_gui(f"NOISeq failed with exit code {proc.returncode}")
            except Exception as e:
                self.show_error_gui(f"Error running NOISeq: {str(e)}")

        threading.Thread(target=worker, daemon=True).start()

    def launch_noisq_annotation(self):
        base_dir = self.params["step1"]["base_output_dir"]
        noisq_results = os.path.join(base_dir, "noisq_results", "noisq_results.xlsx")
        peak_counts = os.path.join(base_dir, "noisq_results", "peak_counts.txt")
        assembly = self.params["step3"]["genome_version"]
        ref_dir = self.params["step3"]["ref_dir"]

        r_script = resource_filename('chromacs', 'annotate_noisq.R')

        cmd = [
            "Rscript", r_script,
            noisq_results,
            peak_counts,
            assembly,
            ref_dir
        ]

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                for line in proc.stdout:
                    self.root.after(0, lambda l=line: self.update_output_gui(l))
                proc.wait()
                self.root.after(0, lambda: messagebox.showinfo(
                    "Success", "NOISeq annotation complete!"
                ))
            except Exception as e:
                self.root.after(0, lambda: self.show_error_gui(str(e)))

        threading.Thread(target=worker, daemon=True).start()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    def run_blocking_command(self, command):
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        # Read output in real-time and update GUI
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                self._update_output(output)

        if process.returncode != 0:
            self.show_error_gui(f"Command failed: {command}")
            return False
        return True

    def update_output_gui(self, text):
        self.root.after(0, self._update_output, text)

    def _update_output(self, text):
        self.output_text.config(state=tk.NORMAL)
        self.output_text.insert(tk.END, text)
        self.output_text.see(tk.END)
        self.output_text.config(state=tk.DISABLED)

    def show_error_gui(self, message):
        self.root.after(0, lambda: messagebox.showerror("Error", message))

    def show_info_gui(self, message):
        self.root.after(0, lambda: messagebox.showinfo("Info", message))


def main():
    root = tk.Tk()
    app = ATACSeqPipeline(root)
    root.mainloop()


if __name__ == "__main__":
    main()