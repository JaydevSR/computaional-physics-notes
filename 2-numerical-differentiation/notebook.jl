### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 7cbdea3c-e478-434d-92d6-4d2ea5627d8b
using PlutoUI

# ╔═╡ 3e749362-258f-4d82-990c-ddc2b9325514
md"# Method of Finite Differences"

# ╔═╡ 3f492f16-f235-4f1b-8802-98bce8b9a9ae
md"
- The  derivative of a function is given by: 

$f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}$

- An approximation will be using a very small value of $h(<< 1)$ for this fraction instead of taking a limit.

- This is an example of the _finite-difference_ approximation and $h$ is called the _step size_.

"

# ╔═╡ 9cae0885-5df6-45cb-9a17-f1bab51c86eb
md"
## Forward Finite Difference

If we expand the function $f(x+h)$ around $x$ with $h$ as the variable we get:

$f(x + h) = f(x) + h f'(x) + \frac{1}{2} h^2f''(x) + \frac{1}{6} h^3 f'''(x) + O(h^4)$

$\implies f'(x) = \frac{f(x+h) - f(x)}{h} + O(h)$

The above is a **first-order** finite difference approximation of the derivative.
"

# ╔═╡ 9484b614-2d27-4b07-8e7c-b459bc2184db
md"
## Backward Finite Difference

If we expand the function $f(x - h)$ instead we get the formula for **backward first-order** finite difference.

$f(x - h) = f(x) - h f'(x) + \frac{1}{2} h^2f''(x) - \frac{1}{6} h^3 f'''(x) + O(h^4)$


$\implies f'(x) = \frac{f(x) - f(x-h)}{h} + O(h)$

- Depends on current function value and a previous function value.
"

# ╔═╡ c52856d7-7cc4-4523-b54d-7ee88990161c
md"
## Central Finite Difference

If we subtract the expansions of $f(x+h)$ and $f(x-h)$ then we get:

$f(x+h) - f(x-h) = 2hf'(x) + \frac{1}{3}h^3f'''(x) + O(h^5)$

$\implies f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)$

- Unlike forward and backward methods, the central finte difference is an second order approximation. Hence, this leads to a closer approximation for same step size.
"

# ╔═╡ bb595b17-e1bd-4a73-a157-210fb580fa58
md"
**_Note on improving the approximations:_** We can obtain higher order approximations from taylor series expansions of $f(x + nh)$ for different integer values of $n$ and then combining those to eliminate powers of $h$. [Table of such approximations](https://en.wikipedia.org/wiki/Finite_difference_coefficient).
"

# ╔═╡ ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
md"
## The Second Derivative

We can also calculate an approximation for the second derivative of a function f(x) using the taylor expansions of $f(x + nh)$.

Using the expansions of $f(x+h)$ and $f(x-h)$ and adding these together we get:

$f (x + h) + f (x − h) = 2f (x) + h^2 f (x) + O(h^4)$

$\implies f''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2} + O(h^2)$

- This is the **second-order central** finite difference approximation to $f''(x)$.
- We can also get other approximations for second derivative using combination of other expansions to improve the precision.
"

# ╔═╡ 8493fa58-4543-49a5-9ecc-0861ce4a3a35
PlutoUI.TableOfContents(title="TOC: Numerical Differentiation")

# ╔═╡ Cell order:
# ╟─7cbdea3c-e478-434d-92d6-4d2ea5627d8b
# ╟─3e749362-258f-4d82-990c-ddc2b9325514
# ╟─3f492f16-f235-4f1b-8802-98bce8b9a9ae
# ╟─9cae0885-5df6-45cb-9a17-f1bab51c86eb
# ╟─9484b614-2d27-4b07-8e7c-b459bc2184db
# ╟─c52856d7-7cc4-4523-b54d-7ee88990161c
# ╟─bb595b17-e1bd-4a73-a157-210fb580fa58
# ╟─ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
# ╟─8493fa58-4543-49a5-9ecc-0861ce4a3a35
