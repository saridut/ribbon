# Spiral equations

## Polar spirals

These spirals are usually described by a function in polar coordinates
$r = f(\theta)$.

### Involute of a circle

_Polar equation._ $r = b \sqrt{1+\theta^2}$, where $b >= 0$ is the radius
of the circle.

_Arc length._ $s = \frac{b}{2} \theta^2 \bigl . \bigr
\rvert_{\theta_0}^{\theta}$, from $\theta_0$ to $\theta$,
 where $\theta > \theta_0$.

_Arc length derivative._ $\frac{ds}{d\theta} = b \theta$

_Curvature._ $\kappa = \frac{1}{b \theta}$. The curvature goes to $\infty$
as $\theta \rightarrow 0$, i.e. at the point where the curve begins on the
circle.

_Tangent._ 
\begin{align}
    t_x &= \frac{\theta \cos \theta - (1+\theta^2) \sin \theta}
                {\sqrt{\theta^2 + (1+\theta^2)^2}}\\
    t_y &= \frac{\theta \sin \theta + (1+\theta^2) \cos \theta}
                {\sqrt{\theta^2 + (1+\theta^2)^2}}
\end{align}

### Generalized Archimedean spiral

_Polar equation._ $r = b \theta^n$, where $b >= 0$ and $n$ is non-zero ($n=0$
is a circle.)

_Arc length._ $s = \frac{b}{n} \theta^n \lvert n\rvert \;
     _2F_1(-\frac{1}{2}, \frac{n}{2}; 1+\frac{n}{2}; -\frac{\theta^2}{n^2})
    \bigl . \bigr \rvert_{\theta_0}^{\theta}$, from $\theta_0$ to $\theta$,
 where $\theta > \theta_0$. 

_Arc length derivative._ $\frac{ds}{d\theta} = 
    b \sqrt{\theta^{2n} + n^2\theta^{2n-2}}$

_Curvature._ $\kappa = \frac{\theta^{1-n}}{b} 
    \frac{\theta^2+n^2+n}{\left(\theta^2+n^2\right)^{3/2}}$. 

_Tangent._ 
```{math}
:label: eq:tgas
\begin{aligned}
    t_x &= \frac{n \cos \theta - \theta \sin \theta}
                {\sqrt{n^2 + \theta^2}}\\
    t_y &= \frac{n \sin \theta + \theta \cos \theta}
                {\sqrt{n^2 + \theta^2}}
\end{aligned}
```

For more details see [Diedrichs (2019)](https://doi.org/10.1017/mag.2019.7).

### Archimedes spiral

_Polar equation._ $r = b \theta$, where $b >= 0$ and $n$ is non-zero ($n=0$
is a circle.)

_Arc length._ From $\theta_0$ to $\theta$, where $\theta > \theta_0$
$
s = \Bigl . \frac{b}{2} \left\{
    \theta \sqrt{1+\theta^2} 
    + \log \left( \theta + \sqrt{1+\theta^2} \right)
    \right \} \Bigr \rvert_{\theta_0}^{\theta}
$.

_Arc length derivative._ $\frac{ds}{d\theta} = 
    b \sqrt{1 + \theta^2}$

_Curvature._ $\kappa = \frac{\theta^2 + 2} 
    {b\left(\theta^2 + 1\right)^{3/2}}$. 

_Tangent._ Set $n = 1$ in equation {eq}`eq:tgas`.

### Fermat spiral

_Polar equation._ $r = b \theta^{1/2}$, where $b >= 0$ and $n$ is non-zero
($n=0$ is a circle.)

_Arc length._ $s = \left . b \sqrt{\theta} \;
    _2F_1\left(-\frac{1}{2}, \frac{1}{4}; \frac{5}{4}; -4\theta^2\right)
      \right \rvert_{\theta_0}^{\theta}$, from $\theta_0$ to
     $\theta$, where $\theta > \theta_0$.

_Arc length derivative._ $\frac{ds}{d\theta} = 
    \frac{b}{2} \sqrt{4 \theta + \frac{1}{\theta}}$

_Curvature._ $\kappa = 2 r \frac{4 r^4 + 3 b^4} 
    {\left(4r^4 + b^4\right)^{3/2}}$. 

_Tangent._ Set $n = \frac{1}{2}$ in equation {eq}`eq:tgas`.

### Hyperbolic spiral

_Polar equation._ $r = b \theta^{-1}$, where $b >= 0$ and $n$ is non-zero
($n=0$ is a circle.)

_Arc length._ $s = \left .
     b \left\{-\frac{\sqrt{1+\theta^2}}{\theta} 
     + \log \left(\theta + \sqrt{1+\theta^2} \right) \right\}
     \right \rvert_{\theta}^{\theta_0}$,
     from $\theta_0$ to $\theta$, where $\theta < \theta_0$.

_Arc length derivative._ $\frac{ds}{d\theta} = 
    b \frac{\sqrt{1+\theta^2}}{\theta^2}$.

_Curvature._ $\kappa = \frac{\theta^4}{b \left(1+\theta^2\right)^{3/2}}$. 

_Tangent._ Set $n = -1$ in equation {eq}`eq:tgas`.

### Lituus

_Polar equation._ $r = b \theta^{-1/2}$, where $b >= 0$ and $n$ is non-zero
($n=0$ is a circle.)

_Arc length._ $s = \left .
    2 \sqrt{\theta}\; _2F_1\left( 
    -\frac{1}{2}, -\frac{1}{4}; \frac{3}{4}; -\frac{1}{4\theta^2} \right)
     \right \rvert_{\theta}^{\theta_0}$,
     from $\theta_0$ to $\theta$, where $\theta < \theta_0$.

_Arc length derivative._ $\frac{ds}{d\theta} = 
    \frac{1}{2\theta} \sqrt{\frac{1}{\theta}+4\theta}$.

_Curvature._ $\kappa = \left(8 \theta^2 - 2\right)
    \left(\frac{\theta}{1+4\theta^2}\right)^{3/2}$. There is a point of
    inflexion at $\theta = \frac{1}{2}$.

_Tangent._ Set $n = -\frac{1}{2}$ in equation {eq}`eq:tgas`.

## Cartesian spirals

These spirals are unsually given in terms of a Césaro equation $\kappa = f(s)$.

### Polynomial spiral

_Césaro equation._ $\kappa = P(s)$, where $P(s)$ is a polynomial in $s$.

_Tangent._ The [Whewell equation](https://en.wikipedia.org/wiki/Whewell_equation)
gives a relation between the tangential angle and the arc length:
 $\frac{d \psi}{ds} = \kappa(s)$, integrating which we obtain 
$\psi = \int_{s_0}^s \kappa(s)$. The tangent given by
\begin{gather} 
t_x = \frac{dx}{ds} = \cos \psi \\
t_y = \frac{dy}{ds} = \sin \psi
\end{gather}

_Cartesian coordinates._ The equations for the cartesian coordinates are
\begin{align}
x &= \int_{s_0}^s \cos \psi\;\mathrm{d}s \\
y &= \int_{s_0}^s \sin \psi\;\mathrm{d}s.
\end{align}

For more details see [Dillen (1990](https://doi.org/10.1007/BF02570761).

### Cornu spiral

_Césaro equation._ $\kappa = a s$.

_Tangent._ 
\begin{align}
    t_x &= \cos \left(\frac{1}{2} a s^2 \right) \\
    t_y &= \sin \left(\frac{1}{2} a s^2 \right)
\end{align}

_Cartesian coordinates._
\begin{align}
    x &= m\; C( s/m ) \\
    y &= m\; S( s/m )
\end{align}
where $m = \sqrt{\pi/a}$ and 
\begin{align}
S(z) &= \int_0^z \sin \left(\frac{\pi t^2}{2}\right)\;\mathrm{d}t\\
C(z) &= \int_0^z \cos \left(\frac{\pi t^2}{2}\right)\;\mathrm{d}t
\end{align}
are the Fresnel integrals. The SciPy implementation is 
<inv:#scipy.special.fresnel>.

### Nielsen spiral

_Césaro equation._ $\kappa = a \exp \left( b s \right)$., where $b$ may be
positive or negative.

_Tangent._ 
\begin{align}
    t_x &= \cos \psi \\
    t_y &= \sin \psi,
\end{align}
where the tangential angle
$\psi = \frac{a}{b} e^{bs} + \frac{\pi}{2}$.
The constant $\frac{\pi}{2}$ comes from integrating the Whewell equation with a
convenient starting point $s_0$.

_Cartesian coordinates._ 
\begin{align}
    x &= -\frac{1}{b}\left(\mathrm{Si} - \mathrm{Si_0} \right) \\
    y &=  \frac{1}{b}\left(\mathrm{Ci} - \mathrm{Ci_0} \right)
\end{align}
where $\mathrm{Si_0}, \mathrm{Ci_0} = \mathtt{sici(a/b)}$ and 
$\mathrm{Si}, \mathrm{Ci} = \mathtt{ sici((a/b)*exp(b*s)) }$.
$\mathtt{sici}$ are `sici` integrals as implemented in the
SciPy function <inv:#scipy.special.sici>.
