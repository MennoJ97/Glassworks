# The held step: literature and theory

Private notes, 2026-11-03. Each claim is marked as either something a source says, with the reference number, or my own inference. Where I could not read a paper I say so. The verification status of every reference is listed at the end.

## 0. The step, stated precisely

The code is `holdingStep` in `web/js/optimizer.js` and its call in `web/js/cavity-design.js`. We minimise $J(x)$, the negated and scaled LDOS, over $x\in[0,1]^n$. The density is $\rho=P_\beta(F(Sx))$: mirror symmetry $S$, cone filter $F$, then second-order subpixel-smoothed projection $P_\beta$. We have $g=\nabla_x J$, $f=\nabla_x k$ and $m=\partial k/\partial\rho$. Our bounds-only MMA step for a gradient $v$ is

$$x_j^+(v)=\mathrm{clip}_{[\alpha_j,\beta_j]}\frac{\sqrt{p_j}\,l_j+\sqrt{q_j}\,u_j}{\sqrt{p_j}+\sqrt{q_j}},\qquad \sqrt{p_j}=(u_j-x_j)\sqrt{v_j^+}+10^{-7},\ \sqrt{q_j}=(x_j-l_j)\sqrt{v_j^-}+10^{-7},$$

with $\alpha_j=\max(0,x_j-\delta,\,l_j+0.1(x_j-l_j))$ and $\beta_j$ defined the same way at the upper end. The held step chooses $\lambda$ so that

$$s(\lambda)\equiv m^{\mathsf T}\big[P_\beta(F(S\,x^+(g+\lambda f)))-\rho\big]=\tfrac12\,(k^\ast-k),$$

where $k$ is the solved pole. $\lambda$ is found by bracketing and then 30 bisections. The next solve is centred at $k+s$. No extra Maxwell solves are needed, only filter and projection passes.

## 1. MMA with constraints (Svanberg 1987, 2002, 2007)

**What the sources say.** The 2007 note [3] defines MMA on problems of the form $\min f_0+a_0z+\sum(c_iy_i+\tfrac12d_iy_i^2)$ subject to $f_i(x)-a_iz-y_i\le 0$ and $x\in X$. Every $f_i$ is replaced by a convex separable approximation $\tilde f_i=\sum_j p_{ij}/(u_j-x_j)+q_{ij}/(x_j-l_j)+r_i$, with

$$p_{ij}=(u_j-x_j)^2\big(1.001(\partial_jf_i)^++0.001(\partial_jf_i)^-+\tfrac{10^{-5}}{x^{\max}_j-x^{\min}_j}\big),$$

and $q_{ij}$ mirrors it. The bounds are $\alpha_j=\max\{x^{\min}_j,\,l_j+0.1(x_j-l_j),\,x_j-0.5(x^{\max}_j-x^{\min}_j)\}$. The asymptotes start at $\pm0.5$ of the range and are then scaled by 0.7 or 1.2 according to whether the variable oscillates, clamped to between 0.01 and 10 times the range [3, eqs. 3.3–3.14]. The subproblem's Lagrangian collects the approximations into one separable function with

$$p_j(\lambda)=p_{0j}+\textstyle\sum_i\lambda_ip_{ij},\qquad q_j(\lambda)=q_{0j}+\sum_i\lambda_iq_{ij}$$

[3, eqs. 5.4–5.5]. The note says the subproblem can be solved by a "dual approach" or by a primal-dual interior-point method. The 2007 Matlab code uses the second [3, §4 end, §5]. The problem form admits only inequality constraints. The note does not discuss equalities. NLopt's MMA documentation says outright that it supports inequality constraints "(but *not* equality constraints)" [17]. The usual route is to write $h\le0$ and $-h\le0$, with the elastic variables $y_i$ absorbing an infeasible subproblem. That last point is my inference from the problem form, not something [3] states.

GCMMA [2, 3] adds inner iterations. If any $\tilde f_i(\hat x)<f_i(\hat x)$ at the trial point, the subproblem is re-solved with larger $\rho_i$, which makes the approximation more conservative. The gradients are not recomputed, but each inner iteration does need the function values $f_i(\hat x)$ [3, §4]. For us those values would be new Maxwell solves.

I could not read the 1987 paper [1] (Wiley paywall). Everything above comes from [3] and the metadata of [1, 2].

**How the held step relates (inference).** With a single constraint, the dual of the MMA subproblem is a one-dimensional concave maximisation over $\lambda$. Its derivative is the approximated constraint value at the primal minimiser $x(\lambda)$. Bisecting on $\lambda$ until the constraint holds is exactly that dual solve. The held step is therefore MMA's dual with $m=1$ and an equality constraint, but it differs in three ways:

1. **The constraint model.** MMA would use the convex rational approximation $\tilde h$ of $k(x)$. We use a model that is linear in $\rho$ and goes through the exact nonlinear filter and projection. Both are first-order accurate at $x$. Ours is more faithful to how the projection bends the step. MMA's version is the convex one.
2. **How the Lagrangian is mixed.** This is the important difference. MMA mixes *coefficients*, $p_j(\lambda)=p_{0j}+\lambda p_{1j}$ (for an equality, use $\lambda<0$ by taking $p,q$ from $-h$). We mix *gradients*: we form $g+\lambda f$ and split that into its positive and negative parts. Where $g_j$ and $\lambda f_j$ have opposite signs, MMA's $p_j(\lambda)$ and $q_j(\lambda)$ are both positive. Its minimiser is then interior and moves continuously with $\lambda$. Ours has one of the two near zero, so it jumps from one bound to the other. **The staircase is an artefact of mixing gradients rather than coefficients. Svanberg's form does not have it.** Even with our stripped-down $p,q$ (no 0.001 or $10^{-5}$ terms), mixing coefficients gives $\sqrt{p_j/q_j}=\sqrt{g_j^+/(\lambda f_j^-)}$ and a smooth $x_j(\lambda)$.
3. **The sign-step behaviour has a cause.** For one variable with a gradient of one sign and no constraint, $\tilde f$ is monotone on $[\alpha_j,\beta_j]$. Its minimiser is then the bound, whatever $|g_j|$ is, so the unconstrained step is $\mathrm{sign}(g_j)$ times a per-variable length set by the asymptotes. Svanberg's 0.001 cross-term does not change this. With symmetric asymptotes at distance $s$, $\sqrt{p/q}=\sqrt{1001}\approx31.6$ gives an interior minimiser at $x-0.939s$, which lies beyond the albefa bound $x-0.9s$. Only the $10^{-5}$ ("raa0") term brings back a dependence on $|g_j|$, and only for $|g_j|\lesssim10^{-2}$ in Svanberg's scaling. Our $10^{-7}$ in square-root space removes even that. So the observation that 81–90 % of variables take the full move, with a cosine of 0.24 to the asked direction, is expected from bounds-only MMA rather than a bug. (The per-variable adaptive length driven by sign changes resembles Rprop [18]. That is an analogy only.)

## 2. Optimality criteria with a bisected multiplier

**What the sources say.** In the 88-line code [5] and its Heaviside companion (the DTU "110 line" file, Nov 2010, read directly), the OC update is `xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))))`. The multiplier `lmid` is found by bisection on $[0,10^9]$ until `(l2-l1)/(l1+l2) <= 1e-3`. **Inside the bisection loop**, the Heaviside version recomputes `xTilde = (H*xnew)./Hs` and `xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta)`, then tests `sum(xPhys(:)) > volfrac*nelx*nely`. Ferrari & Sigmund [7] explain why: the Heaviside projection "is not volume–preserving for any η". It therefore "would require one filter–and–projection application at each bisection step". They also discuss how to cut the number of bisections. The 99-line code [4] has the same multiplier bisection without projection (from the code's description; I did not re-read the 2001 paper itself).

**Relation (inference).** Structurally this is the closest published match to the held step. Both use a closed-form, separable, move-limited update parameterised by one multiplier. Both bisect that multiplier until a constraint, evaluated after the nonlinear filter and projection, takes a set value. Neither needs a new physics solve. The differences:

- the constraint is the frequency, linearised in $\rho$ with $m=\partial k/\partial\rho$, instead of the exactly linear volume;
- the target is half the remaining gap rather than the full constraint;
- the update is MMA's rather than OC's fixed-point heuristic $x_eB_e^{1/2}$;
- our multiplier can have either sign (an equality), whereas OC's is $>0$ (an active inequality).

The volume in OC is held exactly at every iterate. Our frequency is held only to first order, and the half-gap target is what corrects for that.

## 3. Why the Rosen-style projection failed

**Sources.** Rosen [8, 9] projects the gradient onto the null space of the active constraints. For linear constraints a step along that projection stays feasible. For nonlinear constraints Part II adds a move back onto the constraint surface. I could not read either paper (SIAM paywall, no open copy found). This description is the standard one in secondary sources and I have not checked it against the originals.

**Inference.** Rosen's guarantee needs the *step* to lie in $\ker f^{\mathsf T}$. The old code projected the *direction*, $d=g-\frac{g\cdot f}{f\cdot f}f$, and then passed it through a nonlinear map $M$ (the MMA step, then the clamps, then the filter and projection). In general $f^{\mathsf T}M(d)\neq0$. For a near-sign map, $f^{\mathsf T}M(d)\approx-\sum_j\delta_j f_j\,\mathrm{sign}(d_j)$. Nothing makes that vanish: $\sum_j f_jd_j=0$ constrains a weighted sum, not the signs. With a cosine of 0.24 between $M(d)$ and $d$, most of the step lies outside the direction that was asked for, and $f$ sees all of it. The clamped correction $c f$ could not make up the resulting 5–6 bandwidths. The general rule, which is what MMA's dual and OC's bisection both do, is to put the constraint *inside* the subproblem that produces the step, so that it holds for the step actually taken.

## 4. "Half the remaining gap": relaxed linearisation and range-space rate

**Sources.**

- Vardi [10] relaxed the linearised equality to $\gamma_kc_k+G_k\Delta x=0$ with $\gamma_k\in(0,1]$, so that it stays compatible with a trust region. That description comes from Fang et al. [12]. Vardi's own abstract does not mention it, and I could not read the paper.
- Byrd, Schnabel & Shultz [11] minimise a model "subject to a possibly relaxed linearization of the problem constraints" (abstract). Fang et al. [12] say they split the step into a normal part $-\gamma_kG_k^{\mathsf T}(G_kG_k^{\mathsf T})^{-1}c_k$ and a tangential part.
- Feppon, Allaire & Dapogny [13] use the flow $\dot x=-\alpha_J\xi_J-\alpha_C\xi_C$, where the null-space part is $\xi_J=(I-Dg^{\mathsf T}(DgDg^{\mathsf T})^{-1}Dg)\nabla J$ and the range-space (Gauss–Newton) part is $\xi_C=Dg^{\mathsf T}(DgDg^{\mathsf T})^{-1}g$. Along the continuous flow, $g(x(t))=e^{-\alpha_Ct}g(x_0)$ (their eq. 2.16). In the discrete scheme, $g(x_{n+1})=(1-\alpha_C\Delta t)g(x_n)+o(\Delta t)$, and stability requires $0<\alpha_C\Delta t<2$ (App. A, eqs. A.4, A.9). In practice they set $\Delta t=1$ and cap $\alpha_{C,n}\le0.9$, with the range step also limited to a multiple of the mesh size (§5.3.2, eq. 5.14). Bound constraints are handled as ordinary inequality constraints.

**Relation (inference).** Our target $s=\tfrac12(k^\ast-k)$ is Vardi's relaxation with $\gamma=\tfrac12$. It is also Feppon et al.'s geometric decrease with $\alpha_C\Delta t=\tfrac12$, which sits well inside their stable range. Their typical setting is $\le0.9$. The difference is how the two parts are combined. Feppon et al. add a null-space step and a range-space step *linearly*, before any bound handling. The held step instead uses one multiplier on the gradient inside the bound-respecting MMA map, so the range-space correction and the bound handling are solved together. That is the property the Rosen-style version lacked.

## 5. General frame: SLP and SQP with one equality

Inference, with [14] as the general text (not re-read for this note). Because our per-variable model is monotone, $x^+(g+\lambda f)$ minimises $(g+\lambda f)^{\mathsf T}\Delta x$ over the MMA box. The held step is therefore the Lagrangian or dual solution of a **sequential linear program with move limits** and one equality row, $f$-weighted through the projection, with relaxed right-hand side $\tfrac12(k^\ast-k)$. A basic optimal solution of such an LP has at most one variable strictly between its bounds. That single fractional variable is exactly what our bisection cannot represent, which is the staircase seen from the LP side. Solving the LP exactly (sort by $g_j/f_j$, fill in order, leave one variable fractional) would hold the linearised shift exactly. SQP or trust-region SQP would add curvature in place of MMA's asymptotes. Vardi and Byrd–Omojokun-type normal steps are where the "fraction of the gap" idea formally lives.

## 6. Keeping a resonance on target: what others do

- **Liang & Johnson 2013 [15].** They maximise a frequency-averaged LDOS. I could not obtain the full text: Optica needs JavaScript, and the MIT DSpace copy returned HTTP 429. Their abstract says the formulation avoids eigenproblem-based difficulties and allows full 3D topology optimisation. Shaker et al. [16] describe it as maximising LDOS at a fixed $\omega_0$, "perhaps with bandwidth regularization via added loss". As far as I can verify, there is no explicit frequency constraint. The averaging window and the fixed $\omega_0$ pull the mode towards the target. (From memory, unverified here: a Lorentzian average equals the LDOS at a complex frequency $\omega_0+i\Gamma$.)
- **Wang, Christiansen, Yu, Mørk & Sigmund 2018 [19].** They use frequency-domain FEM, MMA and adjoint gradients, maximising a frequency-averaged LDOS at a target wavelength of 1550 nm. The PDF text extracted badly (two overlapping columns), but it mentions "a continuation scheme to attain the final design for the targeted frequency". I could not see any explicit frequency constraint.
- **Shaker, Martinez de Aguirre Jokisch, Chao & Johnson 2025 [16].** This is the closest in spirit. They show that the log-LDOS Hessian has an $O(Q^2)$ eigenvalue whose eigenvector is nearly parallel to $\nabla\mathrm{Re}\,\omega_\ast$, so frequency shifts are what make high-$Q$ LDOS optimisation ill-conditioned. Their fix is to evaluate LDOS at $\mathrm{Re}\,\omega_\ast(p)$, the nearest eigenfrequency found by shift-invert, roughly doubling the number of factorisations. They add a bandwidth *inequality*, often just $\mathrm{Re}\,\omega_\ast\ge\omega_0$, and say it is "not usually active". They optimise with CCSAQ in NLopt, after an unshifted warm-up of up to 1000 iterations. They use SSP [20] in 2D. Compared with us: they let the evaluation frequency follow the mode with an exact eigensolve, while we hold the mode on target with an equality and follow it with a first-order prediction. Both remove the stiff $\nabla k$ direction from what the objective sees. Our "solve at the predicted resonance" is a solve-free approximation of their shift.
- **Structural eigenfrequencies.** There is a long line of work on maximising simple or multiple eigenfrequencies and gaps (Pedersen 2000 [21]; Du & Olhoff 2007 [22]) and on designing for *prescribed* eigenfrequencies (Tsai & Cheng 2013 [23]). I verified only the metadata of these three, not their methods. Ferrari & Sigmund [7] and the OC literature use bisection for *volume*. I found no paper, structural or photonic, that bisects a multiplier inside a closed-form update to put a *frequency* at a set fraction of its gap. One extended search for exactly this came up empty, which is weak evidence that nobody has done it.

## 7. The staircase, quantified

Inference. $s(\lambda)$ has a breakpoint at $\lambda_j=-g_j/f_j$ for every variable. At each breakpoint, variable $j$ jumps by about $\beta_j-\alpha_j\le2\delta$. Through the filter that jump spreads over a disc of radius $r$, so the residual error of the hold is at most about $\max_j|f_j|\cdot2\delta$, in $k$ units. That is small if $\delta\,\|f\|_\infty$ is small compared with the linewidth $k/2Q$. It can be checked directly from the current $f$ and $\delta=\min(0.1,0.2/\beta)$. Bisection is valid wherever $s$ changes sign across the bracket. $s$ is monotone in the linear-in-$x$ model. With SSP it need not be, because the projection depends on $|\nabla\tilde\rho|$ as well as $\tilde\rho$, so a pixel's $\rho$ is not monotone in its neighbours' $x$. The code's fallback (taking the closer end) covers a bracket that fails. Two fixes follow from §1 and §5:

- mix $p,q$ coefficients rather than gradients, which makes $x(\lambda)$ continuous; or
- after bisection, set the one breakpoint variable to its fractional value, which is the LP solution.

## 8. What's new here, if anything

Not much, and that is fine:

- Bisecting a single multiplier inside a separable, closed-form, move-limited update so that a constraint, evaluated through the nonlinear filter and projection, takes a set value is textbook OC as implemented in the DTU codes [5, 7]. It is also MMA's dual with $m=1$ [3].
- Asking for only a fraction of the linearised constraint residual is Vardi's relaxation [10, 11], or a discrete range-space flow with $\alpha_C\Delta t=\tfrac12$ [13].
- Centring the next frequency-domain solve on the moving resonance is the idea of Shaker et al. [16]. They do it with an eigensolve; we do it with a first-order prediction.

The specific combination appears to be unpublished as a named method: a *frequency equality* (not volume), relaxed to half the gap, inside bounds-only MMA, propagated through SSP2, with zero extra solves and the next solve placed at the predicted pole. It is a small engineering variant, not a new algorithm.

It is also slightly *less* principled than plain MMA with one constraint, because it mixes gradients rather than $p,q$ coefficients (§1) and so inherits the staircase. The honest one-line description is: "OC-style multiplier bisection on a linearised frequency equality, Vardi-relaxed by ½, using MMA's separable update."

The one finding that may be worth writing up is negative and pedagogical. Projecting the direction (Rosen) before a sign-like bounded step map does not preserve the constraint. Bounds-only MMA is such a map by construction (§1, item 3).

## References and verification status

DOIs and metadata were checked against Crossref and arXiv IDs against the arXiv API on 2026-11-03. "Read" means I read the full text or the code; "abstract" or "metadata" means only that.

1. K. Svanberg, "The method of moving asymptotes—a new method for structural optimization," *Int. J. Numer. Meth. Eng.* 24(2), 359–373 (1987). doi:10.1002/nme.1620240207. *Metadata only.*
2. K. Svanberg, "A class of globally convergent optimization methods based on conservative convex separable approximations," *SIAM J. Optim.* 12(2), 555–573 (2002). doi:10.1137/S1052623499362822. *Metadata only; GCMMA as described in [3].*
3. K. Svanberg, "MMA and GCMMA – two methods for nonlinear optimization," KTH technical note. The PDF is undated; it describes the "2007 implementations". Formerly at people.kth.se/~krille/mmagcmma.pdf, which is now offline; read via the Internet Archive copy. No DOI. *Read.*
4. O. Sigmund, "A 99 line topology optimization code written in Matlab," *Struct. Multidisc. Optim.* 21(2), 120–127 (2001). doi:10.1007/s001580050176. *Metadata.*
5. E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov, O. Sigmund, "Efficient topology optimization in MATLAB using 88 lines of code," *Struct. Multidisc. Optim.* 43(1), 1–16 (2011; online 2010). doi:10.1007/s00158-010-0594-7. *Code read (top88.m and the Heaviside "110 line" companion from topopt.mek.dtu.dk); paper not read.*
6. M. P. Bendsøe, O. Sigmund, *Topology Optimization: Theory, Methods and Applications*, Springer (2003; Crossref lists 2004). doi:10.1007/978-3-662-05086-6. *Not read for this note.*
7. F. Ferrari, O. Sigmund, "A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D," *Struct. Multidisc. Optim.* 62(4), 2211–2228 (2020). doi:10.1007/s00158-020-02629-w; arXiv:2005.05436. *Read (arXiv).*
8. J. B. Rosen, "The gradient projection method for nonlinear programming. Part I. Linear constraints," *J. SIAM* 8(1), 181–217 (1960). doi:10.1137/0108011. *Metadata only.*
9. J. B. Rosen, "The gradient projection method for nonlinear programming. Part II. Nonlinear constraints," *J. SIAM* 9(4), 514–532 (1961). doi:10.1137/0109044. *Metadata only.*
10. A. Vardi, "A trust region algorithm for equality constrained minimization: convergence properties and implementation," *SIAM J. Numer. Anal.* 22(3), 575–591 (1985). doi:10.1137/0722035. *Abstract; the relaxation is as described in [12].*
11. R. H. Byrd, R. B. Schnabel, G. A. Shultz, "A trust region algorithm for nonlinearly constrained optimization," *SIAM J. Numer. Anal.* 24(5), 1152–1170 (1987). doi:10.1137/0724076. *Abstract.*
12. Y. Fang, S. Na, M. W. Mahoney, M. Kolar, "Fully stochastic trust-region sequential quadratic programming for equality-constrained optimization problems," *SIAM J. Optim.* 34(2), 2007–2037 (2024). doi:10.1137/22M1537862; arXiv:2211.15943. *Read §2 (arXiv), used as a secondary source.*
13. F. Feppon, G. Allaire, C. Dapogny, "Null space gradient flows for constrained optimization with applications to shape optimization," *ESAIM: COCV* 26, 90 (2020). doi:10.1051/cocv/2020015; HAL hal-01972915. No arXiv version found. *Read (HAL v4).*
14. J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer (2006). doi:10.1007/978-0-387-40065-5. *Not re-read.*
15. X. Liang, S. G. Johnson, "Formulation for scalable optimization of microcavities via the frequency-averaged local density of states," *Opt. Express* 21(25), 30812 (2013). doi:10.1364/OE.21.030812. *Abstract only; full text not obtainable.*
16. G. Shaker, B. Martinez de Aguirre Jokisch, P. Chao, S. G. Johnson, "Eigenvalue-accelerated LDOS optimization of high-Q optical resonances," arXiv:2511.16643 (2025). *Read (v1).*
17. NLopt documentation, "NLopt Algorithms," MMA/CCSA section (nlopt.readthedocs.io). *Read.*
18. M. Riedmiller, H. Braun, "A direct adaptive method for faster backpropagation learning: the RPROP algorithm," *Proc. IEEE Int. Conf. Neural Networks*, 586–591 (1993). doi:10.1109/ICNN.1993.298623. *Metadata; cited for the analogy only.*
19. F. Wang, R. E. Christiansen, Y. Yu, J. Mørk, O. Sigmund, "Maximizing the quality factor to mode volume ratio for ultra-small photonic crystal cavities," *Appl. Phys. Lett.* 113(24), 241101 (2018). doi:10.1063/1.5064468; arXiv:1810.02417. *Read (arXiv); the text extraction was partly garbled.*
20. A. M. Hammond, A. Oskooi, I. M. Hammond, M. Chen, S. E. Ralph, S. G. Johnson, "Unifying and accelerating level-set and density-based topology optimization by subpixel-smoothed projection," *Opt. Express* 33(16), 33620 (2025). doi:10.1364/OE.563512; arXiv:2503.20189. *Metadata and abstract.*
21. N. L. Pedersen, "Maximization of eigenvalues using topology optimization," *Struct. Multidisc. Optim.* 20(1), 2–11 (2000). doi:10.1007/s001580050130. *Metadata only.*
22. J. Du, N. Olhoff, "Topological design of freely vibrating continuum structures for maximum values of simple and multiple eigenfrequencies and frequency gaps," *Struct. Multidisc. Optim.* 34(2), 91–110 (2007). doi:10.1007/s00158-007-0101-y. *Metadata only.*
23. T. D. Tsai, C. C. Cheng, "Structural design for desired eigenfrequencies and mode shapes using topology optimization," *Struct. Multidisc. Optim.* 47(5), 673–686 (2013). doi:10.1007/s00158-012-0840-2. *Metadata only.*
24. G. Romano, R. Arrieta, S. G. Johnson, "Differentiating through binarized topology changes: second-order subpixel-smoothed projection," arXiv:2601.10737 (2026). *Metadata and abstract; this is the SSP2 form we use.*

Not verified or not readable: the full texts of [1, 2, 4, 6, 8–11, 14, 15, 21–23], all behind paywalls or unobtainable. Svanberg's KTH homepage now says it is no longer in use, so [3] came from the Internet Archive. No site asked for contact details. Some APIs suggest adding an email, but none was sent.
