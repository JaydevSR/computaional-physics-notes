### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 7cbdea3c-e478-434d-92d6-4d2ea5627d8b
using PlutoUI

# ╔═╡ 8493fa58-4543-49a5-9ecc-0861ce4a3a35
PlutoUI.TableOfContents(title="TOC: Numerical Differentiation")

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

# ╔═╡ 36cbb557-abc3-45f5-8716-930eaf282e25


# ╔═╡ 43f6eb4b-23a4-4f1f-80da-d77ee2fdf458
md"# The Euler Methods"

# ╔═╡ 94c1a797-67a4-4ad6-8a54-10cc7071f5d0


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╟─7cbdea3c-e478-434d-92d6-4d2ea5627d8b
# ╟─8493fa58-4543-49a5-9ecc-0861ce4a3a35
# ╟─3e749362-258f-4d82-990c-ddc2b9325514
# ╟─3f492f16-f235-4f1b-8802-98bce8b9a9ae
# ╟─9cae0885-5df6-45cb-9a17-f1bab51c86eb
# ╟─9484b614-2d27-4b07-8e7c-b459bc2184db
# ╟─c52856d7-7cc4-4523-b54d-7ee88990161c
# ╟─bb595b17-e1bd-4a73-a157-210fb580fa58
# ╟─ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
# ╟─36cbb557-abc3-45f5-8716-930eaf282e25
# ╟─43f6eb4b-23a4-4f1f-80da-d77ee2fdf458
# ╠═94c1a797-67a4-4ad6-8a54-10cc7071f5d0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
