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


# font setup
def load_custom_font(root):
    font_dir = os.path.join(os.getcwd(), "fonts")
    font_path = os.path.join(font_dir, "Roboto-Regular.ttf")

    if os.path.exists(font_path):
        try:
            # Create a new font object and use it
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
            "step5": {}  # Peak Annotation
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

        # Disabling tabs for steps that are not yet completed (keep or remove it??)

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

        # Creating a frame for the output widget and scrollbar
        output_frame = tk.Frame(self.root)
        output_frame.pack(side="bottom", fill="both", expand=True, padx=10, pady=10)

        # Creating a Text widget for output (shows real-time run)
        self.output_text = tk.Text(output_frame, wrap=tk.WORD, height=20, width=100)
        self.output_text.pack(side="left", fill="both", expand=True)

        # Creating a Scrollbar widget
        scrollbar = tk.Scrollbar(output_frame, command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")

        # Linking the Scrollbar to the Text widget
        self.output_text.config(yscrollcommand=scrollbar.set)
        self.output_text.config(state=tk.DISABLED)  # Make it read-only

    def clear_frame(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()

        # required for diffbind
        self.conditions = {}
        self.replicates = {}
        self.controls = {}
        self.control_samples = []  # To store available control samples


    # =================== UI Setup for Each Step =================== #

# Step 1: Selection of raw data and then performing fastqc and multiqc +++++++++++++++++++++++++++++++++++++++++++++++++

    def setup_step1_ui(self):
        # Output directory widgets
        tk.Label(self.step1_frame, text="Base Output Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.output_dir = tk.Entry(self.step1_frame, width=50)
        self.output_dir.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_output_dir, bg= 'grey').grid(row=0, column=2, padx=10, pady=5)

        instruction_text_2 = " [Assign the Base Output Directory where all the corresponding results will be saved. ]"
        tk.Label(self.step1_frame, text=instruction_text_2, wraplength=600, justify=tk.LEFT, anchor="w", font=(self.roboto_font, 9, 'italic')
        ).grid(row=1, column=0, columnspan=3, padx=10, pady=5, sticky="w")

        tk.Label(self.step1_frame, text="Raw Data Directory:", font=(self.roboto_font, 10, 'bold')).grid(row=2, column=0, sticky="w", padx=10, pady=5)
        self.raw_data_dir = tk.Entry(self.step1_frame, width=50)
        self.raw_data_dir.grid(row=2, column=1, padx=10, pady=5)
        tk.Button(self.step1_frame, text="Browse", command=self.browse_raw_data, bg= 'grey').grid(row=2, column=2, padx=10, pady=5)

        instruction_text_1 = "[ If 'Select specific samples' is unchecked, ALL FASTQ files in the RAW DATA DIRECTORY will be processed. ]"
        tk.Label(self.step1_frame, text=instruction_text_1, wraplength=600, justify=tk.LEFT, anchor="w", font=(self.roboto_font, 9, 'italic')
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

        # Save button
        tk.Button(self.step1_frame, text="Save & Next", command=self.save_step1_next, bg="green").grid(
            row=6, column=0, columnspan=3, pady=10)


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
            # Find all FASTQ files with both extensions
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
            match = re.match(r"^(.+?)_[12]\.(?:fastq|fq)\.gz$", base)
            if match:
                sample_name = match.group(1)
            else:
                sample_name = base.replace("_1.fastq.gz", "") \
                    .replace("_2.fastq.gz", "") \
                    .replace("_1.fq.gz", "") \
                    .replace("_2.fq.gz", "") \
                    .split(".")[0]
            sample_names.add(sample_name)

        self.params["step1"]["auto_sample_names"] = sorted(sample_names)

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
                       variable=self.skip_trimming_var, font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0, columnspan=2, padx=10, pady=5)

        # Instructional text
        instruction_3= (
            "Checking this will skip the tirmming done by Trim Galore.\n"
            "Raw data will be used directly for alignment.\n"
            "Please check the FastQC and MultiQC reports if you are unsure of this step. \n"
        )

        tk.Label(self.step2_frame, text=instruction_3, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=5, column=1, columnspan=3, padx=10, pady=5, sticky="w")

        tk.Button(self.step2_frame, text="Save & Next", command=self.save_step2_next, bg="green").grid(row=7, column=1, pady=5)
        tk.Button(self.step2_frame, text="Back", command=lambda: self.notebook.select(0), bg= 'red').grid(row=7, column=0, pady=5)


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
            "Black Bengal Goat (CVASU_BBG_1.0)": {
                "bowtie_ref": "CVASU_BBG_1.0",
                "macs_size": "2.9e9",
                "ensembl_species": "capra_hircus_blackbengal",
                "ensembl_cap": "Capra_hircus_blackbengal",
                "ensembl_assembly": "CVASU_BBG_1.0"
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
        tk.Label(self.step3_frame, text="Select Organism/Genome:", font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, sticky="w", padx=10, pady=5)
        self.genome_var = tk.StringVar()
        self.genome_menu = ttk.Combobox(self.step3_frame, textvariable=self.genome_var,
                                    values=list(self.genome_map.keys()), state="readonly")
        self.genome_menu.grid(row=0, column=1, padx=10, pady=5)

        # Reference status indicator
        self.ref_status = tk.Label(self.step3_frame, text="", fg="red")
        self.ref_status.grid(row=1, column=1, sticky="w", padx=10, pady=5)

        # Buttons
        tk.Button(self.step3_frame, text="Check Reference", command=self.check_reference, bg="gray").grid(row=2, column=1, pady=5)
        tk.Button(self.step3_frame, text="Save & Next", command=self.save_step3_next, bg="green").grid(row=3, column=1,
                                                                                                   pady=5)
        tk.Button(self.step3_frame, text="Back", command=lambda: self.notebook.select(1), bg="red").grid(row=3, column=0,
                                                                                                     pady=5)

    def check_reference(self):
        selected = self.genome_var.get()
        if not selected:
            messagebox.showerror("Error", "Please select a genome reference")
            return False

        genome_info = self.genome_map[selected]

        # Use the script's directory, not __file__
        script_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
        ref_dir = os.path.join(script_dir, "ref_genome", genome_info["bowtie_ref"])

        bt2_base = os.path.join(ref_dir, genome_info["bowtie_ref"])
        fa_file = os.path.join(ref_dir, f"{genome_info['bowtie_ref']}.fa")
        fa_gz = f"{fa_file}.gz"

        # Check index files
        index_files = [f"{bt2_base}.{ext}" for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]
        index_exists = all(os.path.exists(f) for f in index_files)

        if index_exists:
            self.ref_status.config(text="Reference index found!", fg="green")
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

        ref_dir = os.path.join(os.path.dirname(__file__),
                               "ref_genome",
                               genome_info["bowtie_ref"])
        bt2_base = os.path.join(ref_dir,
                                genome_info["bowtie_ref"])

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
        tk.Label(self.step4_frame, text="Select Peak Caller:", font=(self.roboto_font, 10, 'bold')).grid(row=1, column=0, sticky="w", padx=10, pady=5)
        self.peak_caller_choice = tk.StringVar(value="Genrich")
        tk.Radiobutton(self.step4_frame, text="Genrich", variable=self.peak_caller_choice, value="Genrich").grid(row=2,
                                                                                                                 column=0,
                                                                                                                 sticky="w",
                                                                                                                 padx=10)
        tk.Radiobutton(self.step4_frame, text="MACS3", variable=self.peak_caller_choice, value="MACS3").grid(row=2,
                                                                                                             column=1,
                                                                                                             sticky="w",
                                                                                                             padx=10)
        tk.Button(self.step4_frame, text="Save & Next", command=self.save_step4_next, bg="green").grid(row=3, column=2,
                                                                                                      columnspan=3,
                                                                                                      pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.notebook.select(2), bg='red').grid(row=3, column=0, pady=10)

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
        tk.Label(self.step4_frame, text="Select Peak Type:", font=(self.roboto_font, 10, 'bold')).grid(row=1, column=0, sticky="w", padx=10, pady=5)
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
        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0, sticky="w", padx=10, pady=5)
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
        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0, sticky="w", padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=4, column=1, padx=10, pady=5)
        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=4, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(
            row=4, column=3, padx=10, pady=5)

        # Blacklist files
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=5, column=0, sticky="w", padx=10, pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=5, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=5, column=2,
                                                                                                padx=10, pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=5, column=3,
                                                                                                      padx=10, pady=5)

        # Exclude chromosomes
        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=6, column=0, sticky="w", padx=10, pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (mitochondrial chromosomes)", variable=self.exclude_chr_m).grid(row=6, column=1,
                                                                                                padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=6, column=2,
                                                                                                padx=10, pady=5)
        #chromosome exclusion entry
        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=6, column=3, sticky="w", padx=10, pady=5)
        self.custom_chr_entry = tk.Entry(self.step4_frame, width=30)
        self.custom_chr_entry.grid(row=6, column=4, padx=10, pady=5)

        # expand_cut_sites
        tk.Label(self.step4_frame, text='Expand cut-sites to _ bp:', font=(self.roboto_font, 10, 'bold')).grid(row=7, column=0, padx=10, pady=5, sticky='w')
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
        tk.Label(self.step4_frame, text='Assign Max q-value (0 to 1):', font=(self.roboto_font, 10, 'bold')).grid(row=8, column=0, padx=10, pady=5)
        self.max_q_value_var = tk.StringVar(value="0.05")
        self.max_q_combobox = ttk.Combobox(
            self.step4_frame,
            textvariable=self.max_q_value_var,
            values=["0.05", "1"]
        )
        self.max_q_combobox.grid(row=8, column=1, padx=10, pady=5)


        # Merged output name
        self.merged_output_label = tk.Label(self.step4_frame, text="Merged Output Name:", font=(self.roboto_font, 10, 'bold'))
        self.merged_output_label.grid(row=9, column=0, sticky="w", padx=10, pady=5)
        self.merged_output_name_var = tk.StringVar()
        self.merged_output_entry = tk.Entry(self.step4_frame, textvariable=self.merged_output_name_var, width=30)
        self.merged_output_entry.grid(row=9, column=1, padx=10, pady=5)

        # Save button (assign to self.save_button)
        self.save_button = tk.Button(
            self.step4_frame,
            text="Save and Next",
            command=self.save_step4a_settings,
            bg="green"
        )
        self.save_button.grid(row=12, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='red').grid(row=12,
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
        tk.Label(self.step4_frame, text="Test Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=2, column=0, sticky="w", padx=10, pady=5)
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
        tk.Label(self.step4_frame, text="Control Samples:", font=(self.roboto_font, 10, 'bold')).grid(row=3, column=0, sticky="w", padx=10, pady=5)
        self.control_sample_var = tk.StringVar(value="Select a control sample")
        control_dropdown = tk.OptionMenu(self.step4_frame, self.control_sample_var, *step2_samples)
        control_dropdown.grid(row=3, column=1, padx=10, pady=5)

        self.control_listbox = tk.Listbox(self.step4_frame, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.control_listbox.grid(row=3, column=2, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Control Sample", command=self.add_control_sample, bg='gray').grid(
            row=3, column=3, padx=10, pady=5)

        # Instructional text
        instruction_1 = (
            "Select one by one: click on a sample and then click on \'Add Control Sample\'\n"
            "Repeat this for the next sample; you can ADD MULTIPLE SAMPLES at a time  \n"
        )

        tk.Label(self.step4_frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 9, 'italic')
                 ).grid(row=3, column=4, columnspan=4, padx=10, pady=5, sticky="w")

        # Blacklist file selection
        tk.Label(self.step4_frame, text="Select Blacklist Files:", font=(self.roboto_font, 10, 'bold')).grid(row=4, column=0, sticky="w", padx=10, pady=5)
        self.blacklist_listbox = tk.Listbox(self.step4_frame, selectmode=tk.MULTIPLE, height=5, width=50)
        self.blacklist_listbox.grid(row=4, column=1, padx=10, pady=5)
        tk.Button(self.step4_frame, text="Add Files", command=self.browse_blacklist_files, bg='gray').grid(row=4, column=2,
                                                                                                padx=10, pady=5)
        tk.Button(self.step4_frame, text="Remove Selected", command=self.remove_blacklist_files, bg='gray').grid(row=4, column=3,
                                                                                                      padx=10, pady=5)

        # Options for excluding chromosomes
        tk.Label(self.step4_frame, text="Exclude Chromosomes:", font=(self.roboto_font, 10, 'bold')).grid(row=5, column=0, sticky="w", padx=10, pady=5)
        self.exclude_chr_m = tk.BooleanVar(value=False)
        self.exclude_chr_y = tk.BooleanVar(value=False)
        tk.Checkbutton(self.step4_frame, text="Exclude MT (Mitochondrial chromosomes)", variable=self.exclude_chr_m).grid(row=5, column=1,
                                                                                                padx=10, pady=5)
        tk.Checkbutton(self.step4_frame, text="Exclude Y", variable=self.exclude_chr_y).grid(row=5, column=2,
                                                                                                padx=10, pady=5)
        tk.Label(self.step4_frame, text="Custom Chromosomes for Exclusion").grid(row=5, column=3, sticky="w", padx=10, pady=5)
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
        tk.Button(self.step4_frame, text="Save and Next", command=self.save_step4b_settings, bg="green").grid(
            row=9, column=0, columnspan=3, pady=10)
        tk.Button(self.step4_frame, text="Back", command=lambda: self.navigate_back_to_step4_ui(), bg='red').grid(row=9,
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
        instruction_1= (
            "Peak Files will be called and saved in \'peak_files\' directory within the base output directory assigned by the user in step 1\n"
            "The annotated peaks will be saved in 'annotated_dir' directory within the \'peak_files\' directory"
        )
        tk.Label(frame, text=instruction_1, wraplength=600, justify=tk.LEFT, anchor="w",
                 font=(self.roboto_font, 10, 'bold')).grid(row=0, column=0, columnspan=6, padx=10, pady=5, sticky="w")

        # Buttons
        tk.Button(frame, text="Back", command=lambda: self.notebook.select(3), bg='red').grid(row=3, column=0, pady=10)
        tk.Button(frame, text="Run Pipeline", command=self.run_pipeline, bg="green").grid(row=4, column=0, pady=10)

        # Adds DiffBind button
        self.diffbind_btn = tk.Button(
            self.step5_frame,
            text="Run DiffBind Analysis", bg= 'orange',
            command=self.launch_diffbind_config,
            state=tk.NORMAL  # or keep it DISABLED till pipeline finishes??
        )
        self.diffbind_btn.grid(row=5, column=0, pady=10)

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

        try:

            # Get selected samples from Step 1
            selected_samples = self.params["step1"].get("selected_samples")

            # If specific samples are selected, filter Step 2 samples
            if selected_samples:
                # Extract base sample names from selected files (same logic as Step 1)
                selected_names = set()
                for filename in selected_samples:
                    base = os.path.basename(filename)
                    match = re.match(r"^(.+?)_[12]\.(?:fastq|fq)\.gz$", base)
                    if match:
                        sample_name = match.group(1)
                    else:
                        sample_name = base.replace("_1.fastq.gz", "").replace("_2.fastq.gz", "").split(".")[0]
                    selected_names.add(sample_name)

                # Filter Step 2 samples to only include selected ones
                samples = [s for s in samples if s in selected_names]

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
                cmd = f"fastqc -t 8 {files_to_process} -o {fastqc_raw}"
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

            # Step 3: Trim Galore
            # ====== Corrected Trim Galore Check ====== #
            if not self.params["step2"].get("skip_trimming", False):
                self.update_output_gui("Checking Step 3: Trim Galore...\n")
                for sample in samples:
                    # Directly construct expected output paths
                    output_r1 = os.path.join(trimmed_data, f"{sample}_val_1.fq.gz")
                    output_r2 = os.path.join(trimmed_data, f"{sample}_val_2.fq.gz")

                    if not (os.path.exists(output_r1) and os.path.exists(output_r2)):
                        # Get original files from selected_samples
                        sample_files = [f for f in self.params["step1"]["selected_samples"]
                                        if sample in os.path.basename(f)]
                        r1 = [f for f in sample_files if "_1" in os.path.basename(f)][0]
                        r2 = [f for f in sample_files if "_2" in os.path.basename(f)][0]

                        self.update_output_gui(f"Running Trim Galore for {sample}...\n")
                        cmd = (f"trim_galore --basename {sample} --gzip  --cores 2 --paired "
                               f"{r1} {r2} --output_dir {trimmed_data}")
                        if not self.run_blocking_command(cmd):
                            return
                    else:
                        self.update_output_gui(f"Skipping Trim Galore for {sample} (trimmed files exist)\n")
            else:
                self.update_output_gui("Skipping Step 3: Trim Galore (User chose to use raw data)\n")

            # Step 4: FastQC on trimmed data
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
                cmd = f"fastqc -t 8 {trimmed_data}/*.fq.gz -o {fastqc_trimmed}"
                if not self.run_blocking_command(cmd):
                    return
            else:
                self.update_output_gui("Skipping Step 4: FastQC trimmed data (reports exist)\n")

            # Step 5: MultiQC on trimmed data
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
                            f"bowtie2-build --threads 8 {fa_file} {bt2_base}",
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

            # Step 6: Alignment and Sorting
            self.update_output_gui("Checking Step 6: Alignment...\n")
            for sample in samples:
                coord_sorted_bam = os.path.join(bam_output, f"{sample}.sort.bam")
                name_sorted_bam = os.path.join(bam_output, f"{sample}.name.sorted.bam")

                if not (os.path.exists(coord_sorted_bam) and os.path.exists(coord_sorted_bam + ".bai")):
                    self.update_output_gui(f"Aligning {sample}...\n")

                    # Retrieve the correct reference genome (custom or dropdown)
                    genome_version = self.params["step3"].get("genome_version", "").strip()

                    # Use trimmed data if trimming was done, else use raw data
                    if self.params["step2"].get("skip_trimming", False):
                        r1 = f"{raw_data}/{sample}_1.fastq.gz"
                        r2 = f"{raw_data}/{sample}_2.fastq.gz"
                    else:
                        r1 = f"{trimmed_data}/{sample}_val_1.fq.gz"
                        r2 = f"{trimmed_data}/{sample}_val_2.fq.gz"

                    cmd = (f"bowtie2 -p 8 --very-sensitive -X 2000 -x {bt2_base} "
                           f"-1 {r1} -2 {r2} "
                           f"| samtools view --threads 8 -bS - > {bam_output}/{sample}.bam")
                    if not self.run_blocking_command(cmd):
                        return

                    if not os.path.exists(name_sorted_bam):
                        self.update_output_gui(f"Name sorting {sample}...\n")
                        cmd = f"samtools sort -n --threads 8 -o {name_sorted_bam} {bam_output}/{sample}.bam && rm {bam_output}/{sample}.bam"
                        if not self.run_blocking_command(cmd):
                            return
                    else:
                        self.update_output_gui(f"Skipping name sorting for {sample} (Name-sorted BAM exists)\n")

                    self.update_output_gui(f"Coordinate sorting {sample}...\n")
                    cmd = f"samtools sort --threads 8 -o {coord_sorted_bam} {name_sorted_bam} && samtools index -@ 8 {coord_sorted_bam}"
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
                    cmd = f"bamCoverage -p 8 -of bigwig --normalizeUsing=RPKM -v -b {bam_output}/{sample}.sort.bam -o {coverage_bw}"
                    if not self.run_blocking_command(cmd):
                        return
                else:
                    self.update_output_gui(f"Skipping coverage for {sample} (file exists)\n")

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


            # MACS3
            else:
                bam_dir = bam_output
                treated = self.params["step4"].get("treated_samples", [])
                control = self.params["step4"].get("control_samples", [])
                blacklist_files = self.params["step4"].get("blacklist_files", [])
                exclude_chr = self.params["step4"].get("exclude_chr", [])
                max_q = self.params["step4"].get("max_q_value", 0.05)
                genome_size = self.params["step4"].get("genome", "hs")

                cmd_all = ""
                for s in treated:
                    input_bam = f"{bam_dir}/{s}.sort.bam"
                    filtered_bam = f"{bam_dir}/{s}.filtered.sort.bam"

                    self.update_output_gui(f"Filtering BAM for treated sample {s}...\n")
                    try:
                        self.filter_bam(input_bam, filtered_bam, blacklist_files, exclude_chr)
                    except Exception as e:
                        self.show_error_gui(f"BAM filtering failed for sample {s}: {str(e)}")
                        return

                    if control:
                        for ctrl in control:
                            control_input = f"{bam_dir}/{ctrl}.sort.bam"
                            control_filtered = f"{bam_dir}/{ctrl}.filtered.sort.bam"

                            self.update_output_gui(f"Filtering BAM for control sample {ctrl}...\n")
                            try:
                                self.filter_bam(control_input, control_filtered, blacklist_files, exclude_chr)
                            except Exception as e:
                                self.show_error_gui(f"BAM filtering failed for control sample {ctrl}: {str(e)}")
                                return

                            output_prefix = f"{peak_files}/{s}_{ctrl}.macs3.peak"
                            cmd = (f"macs3 callpeak -f BAMPE --call-summits -t {filtered_bam} -c {control_filtered} "
                                   f"-g {genome_size} -n {output_prefix} -B -q {max_q} && rm {filtered_bam} {control_filtered}")
                            cmd_all += cmd + " && "
                    else:
                        output_prefix = f"{peak_files}/{s}.macs3.peak"
                        cmd = (f"macs3 callpeak -f BAMPE --call-summits -t {filtered_bam} "
                               f"-g {genome_size} -n {output_prefix} -B -q {max_q} && rm {filtered_bam}")
                        cmd_all += cmd + " && "

                cmd_all = cmd_all.rstrip(" && ")
                self.update_output_gui("Running Step 8: MACS3 Peak Calling...\n")
                if not self.run_blocking_command(cmd_all):
                    return

            # Step 9: Peak Annotation
            r_script_path = resource_filename('chromacs',
                                              'annotate_peaks_5.R')
            if not os.path.exists(r_script_path):
                self.show_error_gui(f"R script 'annotate_peaks_5.R' not found at:\n{r_script_path}")
                return

            genome_version = self.params["step3"]["genome_version"]
            ref_dir = self.params["step3"]["ref_dir"]
            cmd_annotation = (
                f"Rscript {r_script_path} "
                f"{peak_files} {genome_version} {ref_dir}"
            )

            self.update_output_gui(f"Running Peak Annotation with {r_script_path}...\n")

            if not self.run_blocking_command(cmd_annotation):
                return

            self.root.after(0, self.diffbind_btn.config, {"state": tk.NORMAL})
            self.update_output_gui("\n✅ Pipeline execution complete!\n")

        except Exception as e:
            self.show_error_gui(f"Pipeline failed: {str(e)}")


    # this is required for MACS3 filtering
    def filter_bam(self, input_bam, output_bam, blacklist_files, exclude_chr):

        cmd_parts = []

        # Blacklist filtering
        if blacklist_files:
            blacklist_str = " ".join(blacklist_files)
            cmd_parts.append(f"bedtools intersect -v -abam {input_bam} -b {blacklist_str}")
        else:
            cmd_parts.append(f"samtools view -h -b {input_bam}")

        # Exclude chromosomes
        if exclude_chr:
            pat = "|".join(re.escape(c) for c in exclude_chr)
            cmd_parts.append(
                f"samtools view -h - | awk '$1 ~ /^@/ || $3 !~ /^({pat})$/' | samtools view -b -"
            )

        cmd = " | ".join(cmd_parts) + f" > {output_bam}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise Exception(f"BAM filtering failed: {str(e)}")

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
        self.peak_listbox = tk.Listbox(self.diffbind_window, selectmode=tk.MULTIPLE, exportselection= False, height=10, width=80)
        self.peak_listbox.grid(row=0, column=1, padx=10, pady=5)
        tk.Button(self.diffbind_window, text="Browse Peaks",
                  command=self.populate_peak_list, bg='gray').grid(row=0, column=2)

        # Configuration grid
        self.config_frame = tk.Frame(self.diffbind_window)
        self.config_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=10, sticky='nw')

        # Run button
        tk.Button(self.diffbind_window, text="Run DiffBind",
                  command=self.validate_and_run_diffbind, bg="green").grid(row=3, column=1, pady=10)

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
        ).grid(row=3, column=2, pady=10)

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
            self.controls[name]   = tk.StringVar(value="")

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

        # Validate each row’s required fields
        for idx in self.peak_listbox.curselection():
            peak = self.peak_listbox.get(idx)
            if not self.conditions[peak].get():
                return self.show_error_gui(f"Condition missing for {peak}")
            if not self.replicates[peak].get():
                return self.show_error_gui(f"Replicate missing for {peak}")
            if not self.controls[peak].get():
                return self.show_error_gui(f"Control missing for {peak}")

        self.run_diffbind_analysis()

    def generate_metadata(self):
        rows = []
        base = self.params["step1"]["base_output_dir"]
        for idx in self.peak_listbox.curselection():
            peak_filename = self.peak_listbox.get(idx)
            # Extract SampleID from peak filename
            sample_id = os.path.basename(peak_filename).split('.')[0]

            # Determine PeakCaller based on filename
            if 'macs3' in peak_filename.lower():
                peak_caller = 'narrow'
            else:  # Genrich
                peak_caller = 'bed'

            # Get user-assigned parameters
            condition = self.conditions[peak_filename].get()
            replicate = self.replicates[peak_filename].get()
            control_id = self.controls[peak_filename].get()

            # Build file paths
            bam_reads = os.path.join(base, "bam_output", f"{sample_id}.sort.bam")
            peak_path = os.path.join(base, "peak_files", peak_filename)
            bam_control = os.path.join(base, "bam_output", f"{control_id}.sort.bam") if control_id else ''

            rows.append({
                "SampleID": sample_id,
                "bamReads": bam_reads,
                "Condition": condition,
                "Replicate": replicate,
                "Peaks": peak_path,
                "PeakCaller": peak_caller,
                "ControlID": control_id,
                "bamControl": bam_control
            })

        return pd.DataFrame(rows)

    def run_diffbind_analysis(self):

        # 1. Generate and save the metadata sheet
        df = self.generate_metadata()
        base_dir = self.params["step1"]["base_output_dir"]
        out_csv = os.path.join(base_dir, "diffbind_samplesheet.csv")
        df.to_csv(out_csv, index=False)
        messagebox.showinfo(
            "Metadata Saved",
            f"Metadata sheet written to:\n{out_csv}"
        )

        # 2. Prepare output directory
        out_dir = os.path.join(base_dir, "diffbind_results")
        os.makedirs(out_dir, exist_ok=True)

        # 3. Build the Rscript command
        r_script_path = resource_filename('chromacs', 'diffbind_3.R')
        cmd = ["Rscript", r_script_path, out_csv, out_dir]

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
                        self.update_output_gui("✅ DiffBind analysis completed successfully!\n")
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
