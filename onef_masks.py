"""Source-mask growth shared by the 1/f pilot and stage2-plus reducer."""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy import ndimage as ndi


@lru_cache(maxsize=32)
def _disk(radius: int) -> np.ndarray:
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return xx * xx + yy * yy <= radius * radius


def double_source_mask(source: np.ndarray, max_radius: int = 32) -> tuple[np.ndarray, dict]:
    """Grow each connected source by its equivalent radius, with a safety cap."""
    labels, count = ndi.label(source)
    areas = np.bincount(labels.ravel())
    boxes = ndi.find_objects(labels)
    expanded = source.copy()
    radii = []
    for index, box in enumerate(boxes, start=1):
        if box is None:
            continue
        radius = min(max_radius, max(1, math.ceil(math.sqrt(areas[index] / math.pi))))
        radii.append(radius)
        y0 = max(0, box[0].start - radius)
        y1 = min(source.shape[0], box[0].stop + radius)
        x0 = max(0, box[1].start - radius)
        x1 = min(source.shape[1], box[1].stop + radius)
        patch = labels[y0:y1, x0:x1] == index
        expanded[y0:y1, x0:x1] |= ndi.binary_dilation(patch, structure=_disk(radius))
    return expanded, {"components": count, "radius_cap": max_radius,
                      "median_radius": float(np.median(radii)) if radii else 0.0,
                      "p99_radius": float(np.percentile(radii, 99)) if radii else 0.0}
