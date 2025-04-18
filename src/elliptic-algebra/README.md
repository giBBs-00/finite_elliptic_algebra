# Elliptic Algebra in Python

A modular math library for finite fields, abelian groups, and elliptic curves over \( \mathbb{F}_p \) and \( \mathbb{F}_{p^n} \). Built entirely from scratch — no cryptography libraries, just math and Python. 

---

# Features

- Finite fields \( \mathbb{F}_p \) and extension fields \( \mathbb{F}_{p^n} \)
- Polynomial arithmetic in finite fields
- Cyclic and abelian group structures
- Full elliptic curve arithmetic (point addition, doubling, scalar multiplication)
- Group structure detection: \( \mathbb{Z}/n\mathbb{Z} \), \( \mathbb{Z}/m\mathbb{Z} \times \mathbb{Z}/n\mathbb{Z} \)
- Torsion points and subgroup enumeration
- j-invariant and singularity checks
- Isomorphism between elliptic curve groups and product abelian groups

---

# Project Structure

elliptic-algebra/ ├── finite_field/ │ ├── field.py # Fp and GF (extension fields) │ ├── cyclic.py # Cyclic group structure │ └── abelian.py # AG: finite abelian groups ├── elliptic_curve/ │ └── curve.py # Elliptic curve logic ├── utils/ │ └── helpers.py # Math utilities: gcd, is_prime, etc. ├── examples/ │ └── demo.py # Quick test/demo script ├── .gitignore ├── requirements.txt └── README.md


---

# Example

Run the included example to see it in action:

```bash
python -m examples.demo
```
Output:

The curve equation
All affine points on the curve
The group structure (e.g. Z/36Z or Z/6Z × Z/6Z)
Torsion points of specified order

# Idea

This project is built around classical algebraic structures:

Finite Fields: Basic number systems with prime characteristic
Group Theory: Modeling cyclic and product abelian groups
Elliptic Curves: Algebraic curves used in cryptography and number theory
j-invariant: Classifies elliptic curves up to isomorphism
Torsion: Subgroup of elements of finite order

# Author

Built with curiosity and chaos by a Python beginner.
Ask questions, contribute, or fork for your crypto/math projects.

# License

MIT — use, remix, or break it however you want.