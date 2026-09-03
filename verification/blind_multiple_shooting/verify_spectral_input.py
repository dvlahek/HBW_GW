from pathlib import Path
import hashlib
import numpy as np

PATH = Path(__file__).with_name("blind_q010_spectral_nodes.npz")
EXPECTED_SHA256 = "04016df710591f61023338b20117b9b3fd046f4fc3ca886060f2b4a951d9114f"
EXPECTED_KEYS = {
    "nu_J", "q_JJ", "q_JX", "n_r", "n_phi", "shell_J",
    "nu_H", "q_HH", "shell_H",
}

if not PATH.exists():
    raise SystemExit(f"Missing {PATH.name}")

sha = hashlib.sha256(PATH.read_bytes()).hexdigest()
if sha != EXPECTED_SHA256:
    raise SystemExit(f"SHA256 mismatch: {sha}")

with np.load(PATH, allow_pickle=False) as data:
    keys = set(data.files)
    if keys != EXPECTED_KEYS:
        raise SystemExit(f"Unexpected keys: {sorted(keys)}")

print("SPECTRAL_INPUT_PASS")
print(f"file: {PATH.name}")
print(f"sha256: {sha}")
