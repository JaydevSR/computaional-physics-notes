### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ a4ab362e-f108-11eb-07e1-f199b08c8137
begin
	using PlutoUI
	PlutoUI.TableOfContents(title="Table of Contents")
end

# ╔═╡ 067aa72b-bebb-4873-86e9-b72e7bcbd97c
md"# Preliminaries: Numerical Integration"

# ╔═╡ d262623e-b7cd-42d3-bcc0-c6a365bd4740
md"
## Trapazoidal Approximation

- Given a function $f(x)$ and we need to calculate its integral over an infinitesimal region of length $\delta x$, then we take $F(x)$ as it's antiderivative and we have:

$\int_{x}^{x + \delta x} f(x) dx = F(x + \delta x) - F(x)$
$\int_{x}^{x + \delta x} f(x) dx = F(x) + \delta x F'(x) + \delta x^2 F''(x) / 2 + \delta x^3F'''(x) / 6 + O(\delta x^4)$
$\hspace{0.5in}= \frac{\delta x}{2} [2f(x) + \delta x f'(x) + \delta x^2 f''(x)/3 + O(\delta x^3)]$

- Now we note that this expansion is similar to the taylor expansion for $f(x) + f(x+h)$ with a difference in the fractional coefficient of the third order term. So, we approximate the integral as:

$\int_{x_0}^{x_0 + \delta x} f(x) dx \approx \frac{\delta x}{2}[f(x_0) + f(x_0+h)]$

- This approximation comes with an inherent error that has to do with the difference in third order terms in the approximation and the expansion of integral that is called the _local trapazoidal error_ given by:

$\epsilon_T = \frac{1}{6} \delta x^3 f''(x) - \frac{1}{4} \delta x^3 f''(x) = -\frac{1}{12} \delta x^3 f''(x)$

- Geometrically the trapazoidal approximation is calculating the area of trapazoid bounded by points $\{(x_0, 0), (x_0, f(x_0)), (x_0 + \delta x, f(x_0 + \delta x)), (x_0 + \delta x, 0)\}$
"

# ╔═╡ 6d5545f2-d62f-4606-96c3-ab286dd5d5d5
md"
### Composite Approximation

- Practically when we want to calculate the integral over region a region $[a, b]$ then we devide the region into $N$ intervals of length $\Delta x$ and apply the trapazoidal approximation on these intervals seperately then take the summation:

$\boxed{\int_a^b f(x) dx \approx \frac{1}{2} \sum_{n=0}^{N-1} [f(x + n\Delta x) + f(x+ (n+1) \Delta x)]}$

- If the value of $f''(x)$ is bounded by $M$ in the region then the global error is bounded by: 
$E_G \leq N \times \frac{1}{12} \Delta x^3 M = \frac{b-a}{12} \Delta x^2 M$

- This method is a **second order** method due to the above bound. This method is inherently the same method as Implicit RK2!
"

# ╔═╡ 59e381b9-b313-422d-851c-c3ba91f886a2
html"<br><br>"

# ╔═╡ a232bf5e-ad4e-49e2-aac0-8be838ae9e17
md"## Midpoint Rule"

# ╔═╡ d6978767-6ed6-417f-9ebb-1b90bc867611
md"
- When we make the taylor expansion of the integral of a function $f(x)$ over symmetric interval $\{x_0 - \delta x, x_0 + \delta x\}$ we get:

$\int_{x_0 - \delta x}^{x_0 + \delta x} f(x) dx = \delta x f(x_0) + \frac{1}{24} \delta x^3 f''(x_0) + \mathcal{O}(\delta x^4)$

- Then we make the approximation with a third order error term called _local midpoint error_:

$\int_{x_0 - \delta x}^{x_0 + \delta x} f(x) dx \approx f(x_0) \delta x$

$\epsilon_M = \frac{1}{24} \delta x^3 f''(x)$

- We note that the error is exactly half of the local error in trapazoidal method and similarly the upper bound of global error is about half of this method.
"

# ╔═╡ 516099be-f118-4112-8ec2-b68f0ccc994f
md"### Composite Approximation"

# ╔═╡ 6214c643-e5df-4e86-bd5b-fcb46ef2cb7b
md"
- The composite approximation over an interval $[a, b]$ divided in $N+1$ points i.e. $N$ equal intervals of length $\delta x$ is given by:

$\boxed{\int_a^b f(x) dx \approx \Delta x \sum_{n=0}^{N} f\bigg(a + (n+\frac{1}{2}) \Delta x\bigg)}$

- This method is superior over trapezoidal rule due to half error as well as only one function call per iteration.

**Discrete Grid**

- Taking $x_n = a + n\Delta x$, we can make simplification in the formula as:

$\boxed{\int_a^b f(x) dx \approx \Delta x \sum_{n=0}^{N} f\bigg(\frac{x_{n} + x_{n+1}}{2}\bigg)}$
"

# ╔═╡ 1b209c29-df97-4a9f-9d2a-1259ef1a3d8a
html"<br><br>"

# ╔═╡ 636b73c8-1285-4c2e-98c2-4af68cf82513
md"## Simpson's Rule"

# ╔═╡ c31fd3a0-9cdd-432c-98e5-76c9e3919736
md"
- We note that: $\epsilon_T = -2\epsilon_M$. Hence, we can cancel out the local error from these terms by combining the two methods.

- We take the interval $\{x_0, x_0 + 2 \delta x\}$ and then calculate the integral as:

$I = \int_{x_0}^{x_0 + 2\delta x} f(x) dx$
$\implies I = \frac{2I + I}{3} \approx \frac{2M + T}{3}$

- Here, $M$ and $T$ are midpoint and trapazoidal approximations respectively. So, now finally we have the following approximation:

$I \approx \frac{1}{3}[2(2\delta x f(x_0 + \delta x)) + \frac{1}{2}(2\delta x f(x_0) + 2\delta x f(x_0 + 2 \delta x)) ]$

$\boxed {I \approx \frac{1}{3}\delta x [f(x_0) + 4f(x_0 + \delta x) + f(x_0 + 2\delta x)]}$

- The local error in the simpson's approximation is of order 5. Hence, this is a fourth order method.

- This rule can also be determined by [quadratic interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial) of points $\{x_0, x_0 + \delta x, x_0 + 2 \delta x\}$.
"

# ╔═╡ 9a5fc4c3-3bb0-4e7d-996e-b2fb3ff6c786
md"
### Composite Approximation

In composite form over an interval $[a, b]$ divided in N intervals with the restriction that N has to be even because in each iteration we need three points at each iteration of the method i.e two intervals. So we have:

$\int_{a}^{b} f(x) dx = \frac{\Delta x}{3}[f(a) + 2\sum_{n=1}^{N/2 - 1}f(a + 2 n \Delta x) + 4\sum_{n=1}^{N/2} f(a + (2n - 1) \Delta x) + f(b)]$
"

# ╔═╡ 3c136f2f-090e-4f89-ba28-57e17e2a5de9
html"<br><br>"

# ╔═╡ cce16849-d6e4-4430-8e2e-98fd74d8ecdd
md"## Newton-Cotes Rules"

# ╔═╡ 810c5cc2-f9d0-4eb9-877f-f643b93dd64a
md"
- All of the methods above are part of a family of methods called _Newton-Cotes Rules_. This family consists of methods derived using polynomial interpolations of end points of $n$ equally spaced and connected intervals of infinitesimal size.

- This family of methods has the following general form:

$\int_{x_0}^{x_0 + n\delta x} f(x) dx = \delta x \sum_{i=0}^{n} a_{i, n} f(x_0 + i\delta x)$
$\sum_{i=0}^{n} a_{i,n} = n-1$

- These coefficients $a_{i,n}$ are called **quadrature weights**.

- One example that use quadratic interpolation of 5 points $x_0, \dots, x_0 + 4\delta x$ is [Boole's rule](https://en.wikipedia.org/wiki/Boole%27s_rule) which has a seventh order local error.
"

# ╔═╡ a8a313bf-382b-4add-b579-1157f84f90ff
html"<br>"

# ╔═╡ 4b3114d5-dd47-4673-9c5a-eb514fd85f5b
md"
## Runge's Phenomena

- This refers to the drawback of all the interpolation methods that use intervals of equal length.

- Due to _Runge's phenomena_, Newton-Cotes Rules get increasingly unstable on increasing the number of points used for interpolation. 

- This is due to the fact that such interpolation results in oscillatory functions that diverge more and more from the original function as we move away from center of the interval.

- Newton-Cotes Rules are fairly stable upto 9 point interpolations. After that the stability requirement is not met i.e. all quadrature weights must be non-negetive for stability.

- To overcome this, intervals of different length are used. [More information](https://en.wikipedia.org/wiki/Runge%27s_phenomenon)
"

# ╔═╡ 1cb763b4-2caf-4562-9179-34a56e6618fa
html"<br><br>"

# ╔═╡ 46784270-2f4b-4137-b9a4-ca5567eaeb5c
md"## Gauss-Legendre Quadrature

- If the function to be itegrated is known i.e can be evaluated at any point in the interval of integration then we can use polynomial fit with non-equal spacing between points.

- In fact we can also cluster the points near the end points of interval to prevent fitting any oscillatory functions and thus avoiding Runge's phenomena.

- **Guassian quadrature** is one such method that has several variations depending on the type of quadrature weights.

- **Gauss-Legendre quadrature** method uses quadrature weights derived from _Legendre polynomials_.
"


# ╔═╡ a6ffb3ca-4497-420e-ae5b-00b0ca033c86
md"
### 2-Point Gauss-Legendre Quadrature

- We start with an arbitrary solution using two unknown evaluation points, $x_{-} = x_0 + c_- \delta x$ and $x_{+} = x_0 + c_+ \delta x$. These points lie in the interval $[x_0, x_0 + \delta x]$.

- Thus:

$\int_{x_0}^{x_0 + \delta x} f(x) dx = \frac{1}{2} \delta x [f(x_-) + f(x_+)]$

- Expanding both sides using Taylor series and comparing the coefficients upto fourth order we get the following relations:

$c_- + c_+ = 1 \hspace{0.2in}\text{and}\hspace{0.2in} c_-^2 + c_+^2 = \frac{2}{3}$

$\implies c_{\pm} = \frac{1}{2}\bigg(1\pm \frac{1}{\sqrt{3}}\bigg)$

- Using these values of coefficients we can calculate the integral of the function. The local error is of the **fifth order**. This is the same order as the equivalent Newton-Cotes method i.e Simpson's rule but requires less function calls compared to it.

#### Composite form

- For finite intervals we apply the _2GL_ method over $N$ divisions of the interval to get:

$\boxed{\int_{a}^{b} f(x) dx \approx \frac{1}{2} \Delta x \sum_{n=0}^{N-1}\Bigg[f\bigg( a + (n+\frac{c_-}{2})\Delta x)\bigg) + f\bigg( a + (n+\frac{c_+}{2})\Delta x)\bigg) \Bigg]}$

- This form has a fourth-order global error bound. Hence, this is a **fourth-order method**.
"

# ╔═╡ 8377b796-55e3-4f6c-9fe4-4da544e99bfe
md"**More Information on Guassian Quadrature**: Ref. _section 6.5, Computational Quantum Mechanics by Joshua Izaac and Jingbo Wang_."

# ╔═╡ c94ab624-6374-4bb7-82c9-2b626f9dc858
html"<br><br>"

# ╔═╡ e44b6a31-b464-4916-9f42-d65b7c71d866
md"## Monte Carlo Methods"

# ╔═╡ 9eab867f-f193-4b79-9965-5fc492e91d46
html"<br><br>"

# ╔═╡ 90b569e1-d71a-48a7-bb50-8591a409198c
md"## Packages"

# ╔═╡ ec578cb3-03dd-4377-a5a3-1da004aabfdf
md"
- [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl): Uses Gauss-Konrad method, a Guassian quadrature method.
- **SciPy** and **NumPy**
- [HCubature.jl](https://github.com/stevengj/HCubature.jl) and [NIntegration.jl](https://github.com/pabloferz/NIntegration.jl): For higher dimensional integration.
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
git-tree-sha1 = "94bf17e83a0e4b20c8d77f6af8ffe8cc3b386c0a"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.1"

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
# ╟─067aa72b-bebb-4873-86e9-b72e7bcbd97c
# ╟─a4ab362e-f108-11eb-07e1-f199b08c8137
# ╟─d262623e-b7cd-42d3-bcc0-c6a365bd4740
# ╟─6d5545f2-d62f-4606-96c3-ab286dd5d5d5
# ╟─59e381b9-b313-422d-851c-c3ba91f886a2
# ╟─a232bf5e-ad4e-49e2-aac0-8be838ae9e17
# ╟─d6978767-6ed6-417f-9ebb-1b90bc867611
# ╟─516099be-f118-4112-8ec2-b68f0ccc994f
# ╟─6214c643-e5df-4e86-bd5b-fcb46ef2cb7b
# ╟─1b209c29-df97-4a9f-9d2a-1259ef1a3d8a
# ╟─636b73c8-1285-4c2e-98c2-4af68cf82513
# ╟─c31fd3a0-9cdd-432c-98e5-76c9e3919736
# ╟─9a5fc4c3-3bb0-4e7d-996e-b2fb3ff6c786
# ╟─3c136f2f-090e-4f89-ba28-57e17e2a5de9
# ╟─cce16849-d6e4-4430-8e2e-98fd74d8ecdd
# ╟─810c5cc2-f9d0-4eb9-877f-f643b93dd64a
# ╟─a8a313bf-382b-4add-b579-1157f84f90ff
# ╟─4b3114d5-dd47-4673-9c5a-eb514fd85f5b
# ╟─1cb763b4-2caf-4562-9179-34a56e6618fa
# ╟─46784270-2f4b-4137-b9a4-ca5567eaeb5c
# ╟─a6ffb3ca-4497-420e-ae5b-00b0ca033c86
# ╟─8377b796-55e3-4f6c-9fe4-4da544e99bfe
# ╟─c94ab624-6374-4bb7-82c9-2b626f9dc858
# ╟─e44b6a31-b464-4916-9f42-d65b7c71d866
# ╟─9eab867f-f193-4b79-9965-5fc492e91d46
# ╟─90b569e1-d71a-48a7-bb50-8591a409198c
# ╟─ec578cb3-03dd-4377-a5a3-1da004aabfdf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
