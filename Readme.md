Experimental implementatino of [Crystals-Dilithium signature scheme](https://tches.iacr.org/index.php/TCHES/article/view/839/791)
Work in progress

## Examples

Demonstrate the basic Dilithium API:

```bash
cargo run -p dilithium-core --example basic
# or via the Makefile target
make examples-basic
```

Demonstrate the threshold signing workflow:

```bash
cargo run -p dilithium-core --example threshold
# or via the Makefile target
make examples-threshold
```

Run both examples together:

```bash
make examples
```
