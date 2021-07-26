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
md"# The Euler Method"

# ╔═╡ 94c1a797-67a4-4ad6-8a54-10cc7071f5d0
md"
## The Forward Euler Method

- Given an initial value problem, we can express the ODE in terms of finite differences as:
$\frac{dy}{dx} = f(y(x), x), \hspace{0.1in} y(x_0) = y_0$
$\implies \frac{y(x+h) - y(x)}{h} \approx f(y(x),x)$
$\implies y(x+h) \approx y(x) + hf(y(x),x).$

- In order to discretize the quantities $(y(x), x)$, we choose a step size $h = \Delta x$ and discretize the x-coordinate.

- If we have $N$ discrete values of $x$-coordinate then we can denote these by: $x_0, x_1, \dots, x_{N-1}$ and corresponding $y(x)$ by: $y_0, y_1, \dots, y_{N-1}$.

- So, we have
$\boxed{y_{n+1} \approx y_n + \Delta f(y_n, x_n)}$

- This is an **explicit** method i.e the next point in series is completely determined by the previous points in series. 
"

# ╔═╡ ac9cc650-c65e-4d62-a7df-bad8966c0bc5
md"
## The Backward Euler Method

- If we use the backward finite difference in place of forward, the we get the backward Euler method.

$y'(x) = f(y(x), x)$
$\implies y(x) \approx y(x-h) + hf(y(x), x)$
$\implies y(x+h) \approx y(x) + hf(y(x+h), x+h) \hspace{0.5in} \text{(After change in variables)}$

- After discretizing: 
$y_{n+1} \approx y_n + \Delta x f(y_{n+1}, x_{n+1})$

- This is an **implicit** method i.e. it involves points both at current and next step (_which are unknown_) in calculation. So, for solving this we either have to make it explicit if possible for the given $f(y(x), x)$, or otherwise we rearrange the equation and find the roots of $y_{n+1}$ numerically at each iteration:

$y_{n+1} - y_n - \Delta x f(y_{n+1}, x_{n+1}) = 0$

**Note:** For sufficiently small $\Delta x$ we can take initial guess of root of the above equation to be $y_{n}$ and then make the iterations to better approximate the root because $y_{n+1}$ and $y_{n}$ should be in vicinity.
"

# ╔═╡ c9ab491c-2655-42da-a03e-246470d06fc7
md"
**IMPORTANT THINGS**

- _Step truncation error_: Error that occurs at each step due to appxoimation of taylor series.

- _Global error_: Error that results in final value due to accumulation of step truncation errors. This is equal to the number of iterations times the step truncation error.

- Both forward and backward euler method have second-order step truncation error whereas first order global error. So, both the algorithms are **first-order**.

- We can note that generally backward euler method is more computationally expensive compared to the forward method, as there is an extra step of calculating root at each iteration.

- Backward euler method does have some advantages compared to forward euler method when we talk about stability of the algorithms i.e. algorithm is stable if it doesn't oscillates and grows to infinity.

- Forward euler method has a much restricted set of values of $\Delta x$ for which it is stable compared to backward euler method which is almost always stable except for some $\Delta x$. This can be seem by performing stability analysis. See _5.4, Computational Quantum Mechanics, Joshua Izaax, Jingbo Wang_.
"

# ╔═╡ 266c961c-66b1-48df-975c-0e14b98f7df8


# ╔═╡ 9a17ced8-4c1a-4743-af51-2f45cc318a1d
md"# The Leap Frog Method"

# ╔═╡ 6684910a-0664-474b-ada1-8dd460cf7b74
md"
This method uses central difference formula to solve first order ODEs.

$y'(x) = f(y(x), x)$
$\implies \frac{y(x+h) - y(x-h)}{2h} \approx f(y(x), x)$
$\implies y(x+h) \approx y(x-h) + 2hf(y(x), x)$

In discrete form:

$\boxed{y_{n+1} = y_{n-1} + 2 \Delta xf(y_n, x_n)}$

- In order to use this method explicitly we need two initial conditions i.e. $y_0$ and $y_1$. We can get $y_1$ from $y_0$ using forward euler method.

- The step truncation method of leap frog method is of order $h^3$, and the global error is of order $h^2$. Hence, the leap-frog method is a **second order** method.
- For stability see section _5.5.2 Computational Quantum Mechanics, Joshua Izaax, Jingbo Wang_.
"

# ╔═╡ 6f857855-8b1c-48cf-a56e-20ef3a26b574


# ╔═╡ b3693dc4-b08b-4429-9741-88ced438f5b7
md"
# Packages 

These are the available packages for numerical differentiation:

- [_FiniteDifferences.jl_](https://github.com/JuliaDiff/FiniteDifferences.jl)

- [_FiniteDiff.jl_](https://github.com/JuliaDiff/FiniteDiff.jl)

- [_FDM_](https://github.com/wesselb/fdm)
"

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
# ╟─3e749362-258f-4d82-990c-ddc2b9325514
# ╟─3f492f16-f235-4f1b-8802-98bce8b9a9ae
# ╟─9cae0885-5df6-45cb-9a17-f1bab51c86eb
# ╟─9484b614-2d27-4b07-8e7c-b459bc2184db
# ╟─c52856d7-7cc4-4523-b54d-7ee88990161c
# ╟─bb595b17-e1bd-4a73-a157-210fb580fa58
# ╟─ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
# ╟─36cbb557-abc3-45f5-8716-930eaf282e25
# ╟─43f6eb4b-23a4-4f1f-80da-d77ee2fdf458
# ╟─94c1a797-67a4-4ad6-8a54-10cc7071f5d0
# ╟─ac9cc650-c65e-4d62-a7df-bad8966c0bc5
# ╟─c9ab491c-2655-42da-a03e-246470d06fc7
# ╟─266c961c-66b1-48df-975c-0e14b98f7df8
# ╟─9a17ced8-4c1a-4743-af51-2f45cc318a1d
# ╟─6684910a-0664-474b-ada1-8dd460cf7b74
# ╟─6f857855-8b1c-48cf-a56e-20ef3a26b574
# ╟─b3693dc4-b08b-4429-9741-88ced438f5b7
# ╟─8493fa58-4543-49a5-9ecc-0861ce4a3a35
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
