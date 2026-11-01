"""Upstream Evolution and local Reasoning bindings load together."""
import pyariadne as a
for name in ("Semantics", "Enclosure", "VectorFieldSimulator", "IteratedMapEvolver", "VectorFieldEvolver", "Storage", "LabelledStorage", "ContinuousReachabilityAnalyser", "SafetyCertificate"):
    assert hasattr(a,name), f"Missing Python binding: {name}"
