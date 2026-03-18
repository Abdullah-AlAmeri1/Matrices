# 🧮 Java Matrix Library
 
A clean Java library for matrix operations built from scratch — no external libraries.  
Made as student to practice linear algebra concepts in code.
 
---
 
## 📌 What It Does
 
| Operation | Method |
|---|---|
| Add / Subtract | `add()`, `subtract()` |
| Multiply | `multiply()` |
| Transpose | `transpose()` |
| Scale | `scale()` |
| RREF | `rref()` |
| REF (Row Echelon) | `ref()` |
| Solve linear system Ax = b | `solve()` |
| LU Decomposition | `decomposeLU()` |
| LDU Decomposition | `decomposeLDU()` |
| Determinant | `determinant()` |
| Inverse | `inverse()` |
| Rank | `rank()` |
| Trace | `trace()` |
 
---
 
## 🚀 Getting Started
 
No setup needed. Just copy `Matrix.java` and `Main.java` into your project and start using it.
 

 
## 🧠 Concepts Used
 
- Gaussian & Gauss–Jordan Elimination
- LU and LDU Factorization
- Row Echelon Form (REF) and Reduced Row Echelon Form (RREF)
- Augmented matrices for solving linear systems
- Matrix inverse via `[A | I]` reduction
 
---
 
## ⚠️ Notes
 
- Matrices are represented as `ArrayList<ArrayList<Double>>`
- Throws `IllegalArgumentException` for invalid dimensions
- Throws `IllegalStateException` for singular matrices or inconsistent systems
- Uses `EPS = 1e-9` for floating-point comparisons
 
---
 
## 📚 Why I Built This
 
This is a university project for my **Linear Algebra** course (Year 2).  
The goal was to implement matrix algorithms by hand to really understand how they work — not just call `numpy`.
 
---
*Made with ☕ and a lot of linear algebra homework*
