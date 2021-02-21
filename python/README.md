### Python interface

We export 1 function: `wmh.hash`, which hashes a multiset.
It takes a numpy array of weights, a numpy array of ids, a signature size,
and a sketching method, which defaults to BagMinHash.

Install with `python3 setup.py install`.


```python
import wmh
weights = np.random.rand(10000)
ids = np.random.choice(len(weights))
bmh = wmh.hash(weights, ids, m=1000, stype="bmh")
pmh = wmh.hash(weights, ids, m=1000, stype="pmh1")
assert len(pmh) == 1000
assert len(bmh) == 1000
```
