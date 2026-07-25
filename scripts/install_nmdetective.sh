#!/usr/bin/env bash
# Clone NMDetective-AI (Veiner et al. 2026; Vejni/NMDetectiveAI, MIT) at a pinned commit and fetch
# its Git LFS model weights, then install the package into the active environment. The heavy CUDA
# stack comes from the `nmdetective` Pixi feature; this adds only the GitHub-only package + weights.
set -euo pipefail

SHA="4f647697f46e69b0f975f9817b0ac6767edfea73"
DEST="vendor/NMDetectiveAI"

if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/Vejni/NMDetectiveAI.git "$DEST"
fi
git -C "$DEST" fetch --depth 1 origin "$SHA"
git -C "$DEST" checkout -q "$SHA"
git -C "$DEST" lfs install --local
git -C "$DEST" lfs pull

# The Orthrus remote code the model loads needs this exact HuggingFace stack; the conda solve pins a
# newer huggingface-hub, so pin the trio here (upstream README does the same) before installing.
python -m pip install -U "transformers==4.50.3" "tokenizers==0.21.1" "huggingface-hub==0.30.1"
python -m pip install --no-deps -e "$DEST"
python -c "import NMD; print('NMD package importable')"
ls -la "$DEST/models/NMDetectiveAI.pt"
