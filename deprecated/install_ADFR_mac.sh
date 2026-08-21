#!/usr/bin/env bash
#
# install_adfr_mac.sh
#
# Usage:
#   ./install_adfr_mac.sh [INSTALL_DIR]
#
# If INSTALL_DIR is not provided, defaults to:  ~/Documents/ADFR
#
# This will:
#   0) Choose/install directory
#   1) Download the latest macOS ADFRsuite tarball
#   2) Strip Apple quarantine flags
#   3) Run the ADFRsuite installer into INSTALL_DIR
#
# Notes:
#   - URL is the current macOS ADFRsuite tarball from Scripps.
#   - Script is non-destructive: it will refuse to install into
#     an existing INSTALL_DIR unless you remove it first.

set -euo pipefail

# 0) Install path
INSTALL_DIR="${1:-$HOME/Documents/ADFR}"

echo "=== ADFRsuite macOS installer ==="
echo "Install destination: $INSTALL_DIR"
echo

if [ -e "$INSTALL_DIR" ]; then
  echo "ERROR: Install directory already exists: $INSTALL_DIR"
  echo "       Please remove it or choose a different path."
  exit 1
fi

# 1) Download the macOS ADFRsuite tarball
#    This URL is documented to serve ADFRsuite_x86_64Darwin_1.0.tar.gz for macOS.
ADFR_URL="https://ccsb.scripps.edu/adfr/download/1033/"
TARBALL="adfrsuite_macos.tar.gz"

echo "Step 1: Downloading ADFRsuite tarball from:"
echo "        $ADFR_URL"
echo

if command -v wget >/dev/null 2>&1; then
  wget -O "$TARBALL" "$ADFR_URL"
elif command -v curl >/dev/null 2>&1; then
  curl -L -o "$TARBALL" "$ADFR_URL"
else
  echo "ERROR: Neither wget nor curl is available. Please install one of them."
  exit 1
fi

echo "Download complete: $TARBALL"
echo

# Determine the top-level directory name inside the tarball
ADFR_DIR="$(tar tzf "$TARBALL" | head -1 | cut -d/ -f1)"

echo "Step 2: Extracting tarball into: $ADFR_DIR"
tar xzf "$TARBALL"
echo "Extraction done."
echo

# 2) Strip Apple quarantine flags (if any) from the extracted folder
echo "Step 3: Stripping Apple quarantine attributes (if present)..."
if command -v xattr >/dev/null 2>&1; then
  xattr -dr com.apple.quarantine "$ADFR_DIR" 2>/dev/null || true
else
  echo "Warning: xattr command not found; skipping quarantine stripping."
fi
echo "Quarantine stripping (if needed) complete."
echo

# 3) Run the ADFRsuite installer
echo "Step 4: Running ADFRsuite install.sh..."
cd "$ADFR_DIR"

# -d: destination folder; -c 0: compile .py to .pyc (not .pyo)
./install.sh -d "$INSTALL_DIR" -c 0

echo
echo "=== ADFRsuite installation complete ==="
echo "Installed into: $INSTALL_DIR"
echo "This is your ADFR dir."
