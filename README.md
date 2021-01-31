Contains re-implementations of data structures from https://dl.acm.org/doi/10.1145/3219819.3220089 and https://arxiv.org/abs/1911.00675, as well as Python bindings.


Installation:
```bash
git clone --recursive https://github.com/dnbaker/wmh
cd wmh/python
python3 setup.py install
```

Use:
```python
import wmh
import numpy as np

weights = np.random.rand(1000)
ids = np.random(100000, size=(1000,))
signature_size = 500
sigs = [wmh.hash(weights, ids, m=signature_size, stype=sketcher) for sketcher in
        ("bmh1", "bmh2", "pmh1", "phm1a")]
```

PMinHash effectively normalizes all counts, while bmh treats the id, w pairs as weighted elements in a weighted set.
