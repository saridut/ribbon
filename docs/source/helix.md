# Notes on the helix curve

**Parametric equations** of a right-handed helix with radius $R$,
pitch $P$, and axis along the positive $z$-direction are

\begin{equation}
x = R \cos t, \quad y = R \sin t, \quad z = \left(\frac{P}{2 \pi}\right) t.
\end{equation}

The point $r = \left( x, y, z\right)$ corresponding to $t = 0$ is $\left(
R, 0, 0\right)$.

**Arclength parametrized** form of the above equations is

\begin{equation}
x(s) = R \cos \left( \frac{s}{c}\right), \quad
y(s) = R \sin \left( \frac{s}{c}\right), \quad
z(s) = \frac{bs}{c},
\end{equation}
where $R > 0$ and $c = \sqrt{R^2+b^2}$. The pitch $P = 2\pi b$. 

The equations for the tangent $\bm{T}$, normal $\bm{N}$, and binormal $\bm{B}$
of the Frenet frame are

\begin{align}
\bm{T}(s) & = \left[ 
                -\frac{R}{c} \sin \left(\frac{s}{c}\right),
                 \frac{R}{c} \cos \left(\frac{s}{c}\right),
                 \frac{b}{c}
              \right] \\
\bm{N}(s) & = \left[ 
                -\cos \left(\frac{s}{c}\right),
                -\sin \left(\frac{s}{c}\right),
                 0
              \right] \\
\bm{B}(s) & = \left[ 
                 \frac{b}{c} \sin \left(\frac{s}{c}\right),
                -\frac{b}{c} \cos \left(\frac{s}{c}\right),
                 \frac{R}{c}
              \right].
\end{align}

The curvature $\kappa(s)$ and torsion $\tau(s)$ are

\begin{align}
\kappa(s) &= \lVert \bm{T}'(s)\rVert = \frac{R}{c^2} = \frac{R}{R^2+b^2} \\
\tau(s) &= \frac{b}{c^2} = \frac{b}{R^2+b^2}.
\end{align}

At $s = 0$, $\bs{\beta}(s) = \left[x(s), y(s), z(s)\right] = \left[R, 0,
0\right]$ and

\begin{align}
\bm{T}(0) &= \left[0, \frac{R}{c}, \frac{b}{c} \right]\\
\bm{N}(0) &= \left[-1, 0, 0 \right]\\
\bm{B}(0) &= \left[0, -\frac{b}{c}, \frac{R}{c} \right]\\
\end{align}

For more details see Barrett O'Neill, _Elementary Differential Geometry_ 2nd
ed. p. 60.

If $l$ is the curvature along $s$ and $m$ the twist along $s$,

\begin{gather}
R = \frac{l}{l^2+m^2} \quad \mathrm{and} \quad b = \frac{m}{l^2+m^2} \\
l = \frac{R}{R^2+b^2} \quad \mathrm{and} \quad m = \frac{b}{R^2+b^2}.
\end{gather}

**Helix angle** $\tan \psi = P/2 \pi R \implies P = 2\pi R \tan \psi$.

**Number of turns** $N = L/\sqrt{4\pi^2R^2+P^2}$, where $L$ is the contour
length of the helix.

**Extent** of the helix along its axis $h = NP = PL/\sqrt{4\pi^2R^2+P^2}$.

**Distance** along the axis $\delta = w/\tan \psi = 2 \pi R w/P$.

_Relations between radius/pitch and curvature/twist_.

\begin{gather}
R = \frac{l}{l^2+m^2} \quad \mathrm{and} \quad P = \frac{2\pi m}{l^2+m^2} \\
l = \frac{4\pi^2 R}{4\pi^2 R^2 + P^2} \quad \mathrm{and} 
    \quad m = \frac{2\pi P}{4 \pi^2 R^2 + P^2}.
\end{gather}

**Perpendicularity**. If two helices on a cylinder are at right angles, $\tan
\psi_1 \tan \psi_2 = -1 \implies P_1P_2 = -1 \implies m_1m_2 = -1$.

**Helix axis from quaternion integration**. 
\begin{equation}
\hat{\bm{a}} = \left(-\frac{\upomega_1}{\omag}, 0, \frac{\upomega_3}{\omag} \right),
\end{equation}

where $\bs{\upomega} = \left( l, 0, -m \right)$.

**Translation surface**. Let $H_1$ be a helix with parameters $l_1$ and $m_1$
and a director frame $\left( \bm{d}_1, \bm{d}_2, \bm{d}_3\right)$. Let $H_2$
another helix with handedness opposite to that of $H_1$ with parameters $l_2$
and $m_2$. Then the curvature $n$ along $\bm{d}_1$ at any point on the midline
of the translation surface formed by translating $H_2$ along $H_1$ is $n =
m^2/l$.
