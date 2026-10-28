"""Fig 8 A10-Q4 ink detection rule (image-only; no outcome data). Functions reused by the cohort run."""
from __future__ import annotations
import numpy as np
from scipy import ndimage as ndi


def rgb2lab(rgb: np.ndarray) -> np.ndarray:
    c = rgb.astype(np.float64) / 255.0
    c = np.where(c > 0.04045, ((c + 0.055) / 1.055) ** 2.4, c / 12.92)
    M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])
    xyz = c @ M.T / np.array([0.95047, 1.0, 1.08883]); f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])], -1)


def ink_mask(img: np.ndarray, min_comp: int = 15) -> np.ndarray:
    """v2 (calibrated on the pilot thumbnails only): dark saturated blue/navy, green/teal, black marker; scan-border bars removed."""
    lab = rgb2lab(img); L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]; chroma = np.hypot(a, b); hue = np.degrees(np.arctan2(b, a)) % 360
    navy = (L < 24) & (b < -15) & (hue > 265) & (hue < 320)
    blue = (a < 8) & (b < -22) & (L < 70)
    green = (hue > 90) & (hue < 235) & (chroma > 8) & (L < 80)
    black = (L < 20) & (chroma < 14)
    m = navy | blue | green | black
    H0, W0 = m.shape; bh, bw = max(2, int(0.012 * H0)), max(2, int(0.012 * W0)); m[:bh] = False; m[-bh:] = False; m[:, :bw] = False; m[:, -bw:] = False   # scan border strips
    lab_, n = ndi.label(m, structure=np.ones((3, 3)))
    if n == 0: return m
    H, W = m.shape; keep = np.zeros(m.shape, bool)
    for i, sl in enumerate(ndi.find_objects(lab_), start=1):
        h = sl[0].stop - sl[0].start; w = sl[1].stop - sl[1].start; area = int((lab_[sl] == i).sum())
        if area < min_comp: continue
        border = sl[0].stop >= H - 2 or sl[0].start <= 1 or sl[1].stop >= W - 2 or sl[1].start <= 1
        if border and max(h, w) / max(1, min(h, w)) > 6: continue   # scanner border bar
        keep[sl] |= lab_[sl] == i
    return keep
