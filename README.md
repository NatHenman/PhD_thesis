\title{PhD Thesis codes}
\author{Nathaniel Henman}

\begin{document}

\maketitle

\section{Pre-impact phase of droplet impact onto a viscoelastic surface.}

Sec. 2.1.2.\\
Code: 'pre\_impact\_viscoelastic\_code.m'

\begin{equation*}
    F_{TT}=\frac{1}{\pi} (PV) \int_{-\infty}^\infty \frac{\tilde{e}_3G_{\zeta}(\zeta,T)+\tilde{e}_5G_{\zeta T}(\zeta,T)}{X-\zeta}\mathrm{d}\zeta,
\end{equation*}
\begin{equation*}
    \left((F-G)^3(\tilde{e}_3G_X+\tilde{e}_5G_{XT})\right)_X=12(F-G)_T,
\end{equation*}
\begin{equation*}
    F\rightarrow \frac{X^2}{2}-T, G\rightarrow0,~~~~~\text{as}~T\rightarrow-\infty, |X|\rightarrow\infty.
\end{equation*}

\section{Pre-impact phase of droplet impact onto a flexible surface}

Sec. 2.1.3.\\
Code: 'pre\_impact\_flexible\_code.m'

\begin{equation*}
    F_{TT}=\frac{1}{\pi} (PV) \int_{-\infty}^\infty \frac{P_\zeta(\zeta,T)}{X-\zeta}\mathrm{d}\zeta,
\end{equation*}
\begin{equation*}
    \left((F-G)^3P_X\right)_X=12(F-G)_T,
\end{equation*}
\begin{equation*}    \tilde{e}_1G_{XXXX}+\tilde{e}_2G_{XX}+\tilde{e}_3G+\tilde{e}_4G_{TT}+\tilde{e}_5G_T=P,
\end{equation*}
\begin{equation*}
    F\rightarrow \frac{X^2}{2}-T, G\rightarrow0,P\rightarrow0~~~~~\text{as}~T\rightarrow-\infty, |X|\rightarrow\infty.
\end{equation*}

\section{Pre-impact phase of droplet impact onto a lubricant-infused surface}

Sec. 2.2.\\
Code: 'pre\_impact\_LIS\_code.m'

\begin{equation*}
    F_{TT}=\frac{1}{\pi} (PV) \int_{-\infty}^\infty \frac{P_\zeta(\zeta,T)}{X-\zeta}\mathrm{d}\zeta,
\end{equation*}
\begin{equation*}
    \left(\frac{F^3(F+4\Lambda)}{F+\Lambda}P_X\right)_X=12F_T,
\end{equation*}
\begin{equation*}
    F\rightarrow \frac{X^2}{2}-T,P\rightarrow0~~~~~\text{as}~T\rightarrow-\infty, |X|\rightarrow\infty.
\end{equation*}

\section{Boundary layer jet on a flat lubricant-infused surface}

Sec. 4.2.2.\\
Code: 'BL\_flat\_LIS\_code.m'

\begin{equation*}
    h^2uu_X+2Xh\tilde{v}u_Y=\frac{2X}{\tilde{Re}}u_{YY},
\end{equation*}
\begin{equation}
    (hu)_X+\tilde{v}_Y=0,
\end{equation}
\begin{equation}
    u=\Lambda u_Y,~~~~~\tilde{v}=0,~~~~~\text{at}~Y=0,
\end{equation}
\begin{equation}
    u_Y=0,~~~~~\tilde{v}=0,~~~~~Y=1,
\end{equation}
\begin{equation}
    h=1,~~~~~u=\frac{3}{2+6\Lambda}(2\Lambda+2Y-Y^2),~~~~~\tilde{v}=0,~~~~~\text{at}~X=0
\end{equation}

\section{Boundary layer jet on a deformable lubricant meniscus.}

Sec. 4.3.2.\\
Code: 'BL\_meniscus\_LIS\_code.m'

\begin{equation*}
    h^2X_xuu_X+h\tilde{v}u_Y=\frac{1}{\varepsilon^2\mathrm{Re}}u_{YY},
\end{equation*}
\begin{equation*}
    X_x(hu)_X+\tilde{v}_Y=0,
\end{equation*}
\begin{equation}
    H'''(x)=-\frac{2\mu\mathrm{Ca}}{\varepsilon^3}\frac{u_s(3A+B)}{(H+1)^2}+\frac{\rho\mathrm{Re}\mathrm{Ca}}{\varepsilon}u_s(u_s(A+B+C))',
\end{equation}
\begin{equation*}
    u=\frac{H+1}{\mu h(3A+2B+C)}u_Y,~~~~~\tilde{v}=0,~~~~~\text{at}~Y=0,
\end{equation*}
\begin{equation*}
    u_Y=0,~~~~~\tilde{v}=0,~~~~~\text{at}~Y=1,
\end{equation*}
\begin{equation*}
    u=\frac{3Y}{2}(2-Y),~~~~~\text{at}~X=0,
\end{equation*}
\begin{equation*}
    h\int_0^1u~\mathrm{d}Y=h_0,~~~~~H(0)=H(1)=\int_0^1H~\mathrm{d}x=0.
\end{equation*}

\section{Droplet deformation: flow in air}

Sec. 5.3.1.\\
Code: 'droplet\_air\_code.m'

\begin{equation*}
    \psi_{\xi\xi}+\psi_{\theta\theta}=-e^{2\xi}\zeta,
\end{equation*}
\begin{equation*}
    \zeta_{\xi\xi}+\zeta_{\theta\theta}=\mathrm{Re}(\psi_\theta\zeta_\xi-\psi_\xi\zeta_\theta+e^{2\xi}\zeta_t)
\end{equation*}
\begin{equation*}
    \psi=\psi_\xi=0,~~~~~\text{at}~\xi=0,
\end{equation*}
\begin{equation*}
    \psi\rightarrow2\sinh\xi\sin\theta,~~~~~\zeta\rightarrow\infty,~~~~~\text{as}~\xi\rightarrow\infty,
\end{equation*}
\begin{equation*}
    \psi=\zeta=0,~~~~~\theta=0,\pi.
\end{equation*}

\section{Droplet deformation: flow in droplet}

Sec. 5.3.2.\\
Code: 'droplet\_inside\_code.m'

\begin{equation*}
    r^2P_{rr}-rP_r+P+P_{\theta\theta}=0
\end{equation*}
\begin{equation*}
    r^2U_t=-r^2P_r+rP+\frac{1}{\mu_1\mathrm{Re}}(r^2U_{rr}-rU_r+U_{\theta\theta}-2V_\theta)
\end{equation*}
\begin{equation*}
    r^2V_t=-rP_\theta+\frac{1}{\mu_1\mathrm{Re}}(r^2V_{rr}-rV_r+V_{\theta\theta}+2U_\theta)
\end{equation*}
\begin{equation*}
    U=V=P=0,~~~~~\text{at}~r=0,
\end{equation*}
\begin{equation*}
    P=p_g|_{r=1}+\frac{1}{\mathrm{We}},~~~~~V_r-V=\mu_1\zeta_g|_{r=1}-U_\theta,~~~~~U_r=-V_\theta,~~~~~\text{at}~r=1,
\end{equation*}
\begin{equation*}
    U_\theta=V=P_\theta=0,~~~~~\text{at}~\theta=0,\pi.
\end{equation*}




