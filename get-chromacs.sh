#!/bin/bash
set -e

ENV_NAME="chromacs"
REPO_URL="https://github.com/epigen-bioinfolab/CHROMACS_V1.git"
INSTALL_DIR="$HOME/.local/bin"
MIN_CONDA_VERSION="23.1"

# ---------------------------
# Check for Conda
# ---------------------------
echo "🔍 Checking for Conda..."
if ! command -v conda &> /dev/null; then
  echo "❌ Conda not found. Please install Miniconda or Mamba first:"
  echo "   https://docs.conda.io/en/latest/miniconda.html"
  exit 1
fi

# ---------------------------
# Version Check
# ---------------------------
CURRENT_VERSION=$(conda --version | awk '{print $2}')
VERSION_OK=$(python3 -c "from packaging import version; print(version.parse('$CURRENT_VERSION') >= version.parse('$MIN_CONDA_VERSION'))" 2>/dev/null || echo "False")

if [ "$VERSION_OK" != "True" ]; then
  echo "⚠️ Your Conda version ($CURRENT_VERSION) is older than required ($MIN_CONDA_VERSION)."
  echo "👉 Please update Conda before proceeding:"
  echo "    conda update -n base -c defaults conda"
  exit 1
fi

# ---------------------------
# Prepare environment
# ---------------------------
echo "📂 Cloning ChroMACS repository..."
git clone "$REPO_URL" "$HOME/CHROMACS_INSTALL_TEMP"
cd "$HOME/CHROMACS_INSTALL_TEMP"

eval "$(conda shell.bash hook)"

echo "🧪 Creating Conda environment: $ENV_NAME"
conda env create -f environment.yml -n "$ENV_NAME" || {
  echo "⚠️ Environment '$ENV_NAME' may already exist. Trying to activate..."
}
conda activate "$ENV_NAME"

echo "📦 Installing ChroMACS with pip..."
pip install .

# ---------------------------
# Create shortcut
# ---------------------------
mkdir -p "$INSTALL_DIR"
echo -e "#!/bin/bash\nconda activate $ENV_NAME && python -m chromacs \"\$@\"" > "$INSTALL_DIR/chromacs"
chmod +x "$INSTALL_DIR/chromacs"

if ! echo "$PATH" | grep -q "$INSTALL_DIR"; then
  echo "⚠️  $INSTALL_DIR is not in your PATH."
  echo "👉 Add this to your ~/.bashrc or ~/.zshrc:"
  echo "export PATH=\"\$PATH:$INSTALL_DIR\""
fi

# ---------------------------
# Done
# ---------------------------
echo -e "\n✅ ChroMACS installed successfully!"
echo "👉 To activate the environment manually: conda activate $ENV_NAME"
echo "👉 To run from anywhere: chromacs [options]"

rm -rf "$HOME/CHROMACS_INSTALL_TEMP"

