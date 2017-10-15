The Physics of CGMF
*******************

:program:`CGMF` follows the de-excitation of primary fission fragments through the evaporation of neutrons and photons until they reach a ground-state or a long-lived excited state. The statistical Hauser-Feshbach theory is used to infer the probabilities of evaporating a neutron or a photon at each step of the decay. Some details on the physics used is described in this Chapter.


Initial Fission Fragment Yields
===============================

:program:`CGMF` does not calculate the initial pre-neutron fission fragment yields. Instead, it reads or reconstructs those yields in mass, charge and kinetic energy, :math:`Y(A,Z,TKE)`, from experimental data or systematics. Several theoretical efforts are underway to predict fission fragment yields from dynamical fission calculations. We will incorporate the results of those works as they become available.

In the present version of :program:`CGMF`, only binary fission events are considered. Ternary fission where an alpha-particle is emitted along with the two fragments is not treated, nor more complicated "fission" splitting, e.g., accompanied with cluster emission. In addition, the neutron emission is assumed to happen only once both fragments are `fully accelerated`. In other words, no `scission` neutrons are considered at this point. However, multi-chance fission processes such as (n,n'f), (n,2nf), etc., as well as pre-equilibrium contributions are taken into account at higher incident energies.


Mass Yields
-----------

Thermal Neutrons and Spontaneous Fission
++++++++++++++++++++++++++++++++++++++++

For important fission reactions such as the thermal neutron-induced fission cross-section of Pu-239 and U-235, enough reliable experimental data exist to reconstruct those initial yields reasonably well. This was done for instance in the case of thermal neutron-induced fission on Pu-239 in [Talou:2011]_. In that case, the fission fragment mass distribution :math:`Y(A)` was obtained from a least-square fit of several experimental data sets, as shown in Fig. MassYields_.


.. _MassYields:

.. figure:: _static/images/Pu239T_YA.pdf
    :width: 85%

.. figure:: _static/images/Pu239T_YZ.pdf
   :width:  85%
	
   Primary fission fragment mass (top) and charge (bottom) yields for thermal neutron-induced fission of Pu-239. Experimental data on the mass yields were used in a least-square fit to produce the black line. The charge distribution was reconstructed following the Wahl systematics for each fragment mass, as explained below.



Incident Neutron Energies up to 20 MeV
++++++++++++++++++++++++++++++++++++++

At higher incident neutron energies, experimental data become scarce or non-existent, and one has to rely on theoretical models to construct the fragment yields. In the version 1.0.6 of the code, we have implemented a simplified energy dependence for the mass yields. It consists in using a three Gaussian model, whose parameters have been adjusted to reproduce experimental data, when available. For a particular incident neutron energy :math:`E_n`, the yield for the fragment mass :math:`A` is given by:

.. math::

	Y(A;E_n) = G_0(A)+G_1(A)+G_2(A),

where :math:`G_0` corresponds to a symmetric mode,

.. math::

	G_0(A)=\frac{W_0}{\sigma_0\sqrt{2\pi}}\mbox{exp}\left(-\frac{(A-\overline{A})^2}{2\sigma_0^2}\right),

and :math:`G_1` and :math:`G_2` to two asymmetric modes

.. math::

	G_{1,2}(A) = \frac{W_{1,2}}{\sigma_{1,2}\sqrt{2\pi}} \left[ \mbox{exp}\left(-\frac{(A-\overline{A}-D_{1,2})^2}{2\sigma_{1,2}^2}\right) + \mbox{exp}\left(-\frac{(A-\overline{A}+D_{1,2})^2}{2\sigma_{1,2}^2}\right) \right].

Here, :math:`\overline{A}=A_f/2` with :math:`A_f` the mass of the fissioning system, which can differ from the original compound nucleus if pre-fission neutrons are emitted. The parameters :math:`D_i` are governed by spherical and deformed shell closures. Their values decrease by 1/2 for each pre-fission neutron emitted. The energy-dependence for the width parameters is given by:

.. math::

	\sigma_i = \sigma_i^{(0)}+\sigma_i^{(1)}E_n+\sigma_i^{(2)}E_n^2

for :math:`i=1,2`. The width of the symmetric mode :math:`\sigma_0` is assumed to be energy independent.

The weights :math:`W_i` of the Gaussians depend slowly on the incident energy, with an increasing symmetric component. For :math:`W_{1,2}`, we adopt the following energy dependence:

.. math::

	W_i = \frac{W_i^0}{1+\mbox{exp}[(E_n-E_1)/E_2]},

with two adjustable parameters :math:`E_{1,2}`. The weight :math:`W_0` for the symmetric mode is obtained through the normalization condition

.. math::

	W_0 + W_1 + W_2 = 2.

.. _fig_YAKE-Einc:

.. figure:: _static/images/YAKE_Einc.pdf
    :width: 60%
    :align: center

    Fission fragment yields as a function of mass and kinetic energy, for several incident neutron energies in the neutron-induced fission reaction on Pu-239. Multi-chance fission and pre-equilibrium contributions are taken into account as the incident neutron energy increases.

If neutrons are emitted prior to fission, the fissioning nucleus is formed with a residual excitation energy smaller than the initial excitation energy. In this case, an "equivalent" incident neutron energy is defined as the neutron energy that would produce the :math:`(A_0-\nu_{pre})` fissioning nucleus, with :math:`\nu_{pre}` pre-fission neutrons, at the same residual excitation energy. Hence, :math:`E_n` becomes

.. math::

	E_n = E^*-S_{n|A_0-\nu_{pre}}.

The same equivalent incident energy is used in the Wahl parameterization for the charge distribution.

In the current version of the code, we impose that :math:`E^*` be greater or equal than the fission barrier height in the :math:`(A_0-\nu_{pre})` nucleus, and therefore neglect any subbarrier fission events.

.. note::

	Initial parameterizations for the three-Gaussian model were taken from the :program:`FREYA` code. Newer parameterizations based on better fits to known experimental data are being investigated.


Charge Yields
-------------

Wahl systematics [Wahl:2002]_ are then used to obtain the charge distribution for a given mass following:

.. math::
   :label: YZA

	P(Z|A) = \frac{1}{2}F(A)N(A)\left[ erf(V)-erf(W) \right],

where

.. math::
	V = \frac{Z-Z_p+0.5}{\sigma_z\sqrt(2)} \mbox{ and } W=\frac{Z-Z_p-0.5}{\sigma_z\sqrt(2)}

and :math:`erf(x)` represents the error function. The factor :math:`N(A)` is simply a normalization factor. The most probable charge is given by

.. math::
   :label: Zp

	Z_p=A_h\frac{Z_c}{A_c}+\Delta Z,

where :math:`Z_c,A_c` are the charge and mass of the fissioning compound nucleus, :math:`\sigma_z` is the charge width parameter and :math:`\Delta Z` is the charge deviation. The odd-even factor :math:`F(A)` is computed as

.. math::
   :nowrap:

	\begin{alignat*}{3}
	F(A) &= F_Z\times F_N && \mbox{for $Z$ even and $N$ even} \nonumber \\
	F(A) &= F_Z/F_N && \mbox{for $Z$ even and $N$ odd} \nonumber \\
	F(A) &= F_N/F_Z && \mbox{for $Z$ odd and $N$ even} \nonumber \\
	F(A) &= 1/(F_Z\times F_N) && \mbox{for $Z$ odd and $N$ odd} \nonumber
	\end{alignat*}

The average charge distribution is obtained by convoluting :math:`Y(Z|A)` over the fragment mass distribution :math:`Y(A)`, and the result is shown in figure fig-YZ-Einc_ for the heavy fission fragments only.

.. _fig-YZ-Einc:

.. figure:: _static/images/YZ_Einc.pdf
   :width: 70%
   :align: center

   Fission fragment charge distribution as a function of incident neutron energy for the Pu-239 (n,f) reaction.


Total Kinetic Energy (TKE) Distributions
----------------------------------------

The average total kinetic energy :math:`\overline{TKE}` is an important quantity that determines in great part the total excitation energy available in the system for the evaporation of neutrons and photons. Since most neutrons are emitted prior to photon emission, the average total prompt neutron multiplicity, :math:`\overline{\nu}`, strongly depends on an accurate value for :math:`\overline{TKE}`. For the simulation of single fission events, :math:`TKE` distributions have to be known for all fragments.

For thermal neutron-induced fission reactions on important isotopes as well as spontaneous fission, some reliable and rather consistent experimental data exist, albeit less so in the symmetric region where fission events are rare.

To reconstruct the total kinetic energy dependence of the fission fragment yields, one can use experimental information on the average :math:`TKE` as a function of the fragment mass :math:`A` as well as its width :math:`\sigma_{TKE}(A)`. Continuing on the example above for thermal neutron-induced fission of Pu-239, we have performed a least-square fit of :math:`\overline{TKE}(A)` as seen in Fig. fig-TKEA_.

.. _fig-TKEA:

.. figure:: _static/images/Pu239T_TKE_A.pdf
   :width:  70%
   :align:  center

   Average total kinetic energy as a function of the heavy fragment mass in the case of the thermal neutron-induced fission of Pu-239.

The :math:`TKE` distribution for each fragment mass is then reconstructed using

.. math::

	P(TKE|A) = \left( 2\pi \sigma^2_{TKE}(A) \right)^{-1/2} \times \exp\left[ -\frac{\left[ TKE-\overline{TKE}(A)\right]^2}{2\sigma^2_{TKE}(A)} \right].

In a first approximation, one can assume that the shape of :math:`\overline{TKE}(A)` as well as :math:`\sigma_{TKE}(A)` are independent of the particular fissioning system and the energy of the incident neutron (see Fig. fig-TKEA-Isotopes_). We therefore assume that only the absolute scaling of :math:`\overline{TKE}` changes with energy.

.. _fig-TKEA-Isotopes:

.. figure:: _static/images/TKEvsA.png
	:width:  60%
	:align:  center

.. figure:: _static/images/sigTKEvsA-U238.png
	:width:  60%
	:align:  center

	Experimental data available for the mass and incident energy dependence of :math:`\overline{TKE}` and :math:`\sigma_{TKE}` are shown for several fissioning systems and incident neutron energies.

.. note::

	The mass-dependent average total kinetic energy does change with incident energy, reflecting changes in the shell corrections as the excitation energy is increased. A more refined treatment of this quantity will be tackled in the future.

The energy-dependence of :math:`\overline{TKE}` is poorly known for most systems. However, recent experimental data have shed some light on this issue. In the current version of the code, we assume that for each pair of fission fragments, :math:`TKE` can be represented by a normal distribution :math:`\mathcal{N}_{(\langle TKE \rangle,\sigma_{TKE})}(A,E_n)`, and assume that the energy dependence is entirely encoded in the average value :math:`\overline{TKE}`. 

In the current code implementation, the mass and energy-dependent distributions :math:`TKE(A,E_n)` are obtained as

.. math::

	\overline{TKE} (A,E_n) = \overline{TKE} (A,E_{th}) \times \frac{ \overline{TKE}(E_n)}{\sum_A{Y(A,E_n)\overline{TKE}(A,E_{th})}}

The energy dependence of :math:`\overline{TKE}(A)` is given by the Madland systematics [Madland:2006]_, which are simple linear or quadratic fits to experimental data for selected isotopes. Making the distinction between the total fission fragment (pre-neutron) kinetic energy, :math:`TKE_{pre}`, and the total fission product (post-neutron) kinetic energy, :math:`TKE_{post}`, those systematics read:

For **n+U-235**,

.. math::
   :nowrap:

	\begin{eqnarray}
	TKE_{pre} &=& (170.93\pm0.07)-(0.1544\pm0.02)E_n \mbox{ (MeV)}, \nonumber \\
	TKE_{post} &=& (169.13\pm0.07)-(0.2660\pm0.02)E_n \mbox{ (MeV)}.
  	\end{eqnarray}

For **n+U-238**,

.. math::
   :nowrap:

	\begin{eqnarray}
	TKE_{pre} &=& (171.70\pm0.05)-(0.2396\pm0.01)E_n + (0.003434\pm0.0004)E_n^2 \mbox{ (MeV)}, \nonumber \\
	TKE_{post} &=& (169.8\pm0.05)-(0.3230\pm0.01)E_n + (0.004206\pm0.0004)E_n^2 \mbox{ (MeV)}.
  	\end{eqnarray}
 
And for **n+Pu-239**,

.. math::
   :nowrap:

	\begin{eqnarray}
	TKE_{pre} &=& (177.80\pm0.03)-(0.3489\pm0.02)E_n \mbox{ (MeV)}, \nonumber \\
	TKE_{post} &=& (175.55\pm0.03)-(0.4566\pm0.02)E_n \mbox{ (MeV)}.
  	\end{eqnarray}

Madland's fits were only constructed up to the threshold for second-chance fission. We assume however that they are valid at higher energies as well for the initial fissioning nucleus. Above the second-chance fission threshold, the average :math:`TKE` does not necessarily follow a linear or quadratic behaviour though, as successive neutron emissions modify the fissioning nucleus and its excitation energy. We further assume that Madland's energy-dependence parameterizations remain valid for the nuclei A-1, A-2, etc. Only the reference thermal value of :math:`\overline{TKE}(E_{th})` is changed according to Viola's systematics [Viola:1985]_

.. math::
  :label: Viola

	\overline{TKE}_{th} = (0.1189\pm0.011)\frac{Z^2}{A^{1/3}}+(7.3\pm1.5) \mbox{ MeV}.


.. _fig_YKE_Einc::

.. figure:: _static/images/YKE_Einc.pdf
   :width: 70%
   :align: center

   Fission fragment kinetic energy distribution as a function of incident neutron energy for the Pu-239 (n,f) reaction.


Complete :math:`Y(A,Z,TKE)` Yields Reconstruction
-------------------------------------------------

Finally, the full pre-neutron emission fission fragment distributions can be reconstructed as:

.. math::
   :label: YAZTKE

	Y(A,Z,TKE) = Y(A) \times P(Z|A) \times P(TKE|A) 

The resulting :math:`Y(A,TKE)` distribution is shown here:

.. _fig-YATKE:

.. figure:: _static/images/Pu239T_YATKE.png
   :width:  70%
   :align:  center

   Mass and Total Kinetic Energy yields reconstructed using Eq. :eq:`YAZTKE` in the thermal neutron-induced fission of Pu-239.

The approach described above to evaluate the pre-neutron emission fission fragment yields is not unique, and depends on the type of experimental data that have been measured. In some cases, the two-dimensional :math:`Y(A,TKE)` distribution has been measured [Hambsch:2007]_ [Romano:2010]_, and therefore only the charge distribution for every fragmentation has to be computed to obtain the full distribution. In the majority of cases, however, no such information is available and one has to rely on systematics and/or phenomenological models. The present version of :program:`CGMF` is limited to the few isotopes and reactions that have been well measured. The extension to other isotopes and reactions is planned for the near future.

Pre-Fission Neutrons
====================

If the initial excitation energy in the compound nucleus is high enough, there is a chance that neutrons are evaporated prior to fission. We then talk about first-chance :math:`(n,f)`, second-chance :math:`(n,n'f)`, third-chance :math:`(n,2nf)`, etc., fissions. The probabilities for each multi-chance fission event to occur can be computed from the :math:`\Gamma_n/\Gamma_f` ratio as a function of the incident neutron energy. This ratio depends in turn on the fission barrier heights in the various compound nuclei :math:`A, A-1, A-2`, etc. The :program:`CoH-3.0.4` code was used to calculate those ratios for different actinides. As an example, we show here the case of n+Pu-239, in comparison with ENDF/B-VII.1 and JENDL-4.0 evaluations. The :program:`CoH` calculations tend to predict a much higher second-chance fission probability at the expense of the first-chance, compared to the evaluations. These quantities are not observables though, and it is therefore difficult to judge about the validity of those curves at this point.

.. figure:: _static/images/Pu239-multichancefission.pdf
  :width: 75%
  :align: center

  Multi-chance fission probabilities in the neutron-induced fission reaction on Pu-239 as calculated with the :program:`CoH` code (and used in :program:`CGMF`), and in comparison with the ENDF/B-VII.1 and JENDL-4.0 evaluations.

In :program:`CGMF`, those multi-chance fission probabilities are sampled to determine the number of pre-fission neutrons. Then, the energies of those neutrons are obtained by sampling the corresponding neutron spectra. In the case of the first emitted neutron, the spectrum corresponds to a weighted sum of a pre-equilibrium and an evaporation components. The fraction of pre-equilibrium neutrons is also calculated in the :program:`CoH` code using the exciton model. Then, the first neutron-out spectrum is given by:

.. math::

  \chi_1 = f_{pe}\chi_{pe}+(1-f_{pe})\chi_{evap}.

The energy-dependent fraction :math:`f_{pe}` can be fitted by a simple function:

.. math::

  f_{pe}(E_{inc}) = \frac{1}{1+\exp\left[ (12.49-E_{inc})/10.21 \right]}-0.042 E_{inc} -0.25.

As can be seen in Fig. fig-PE_, it is a very reasonable approximation for neutron-induced reactions on U-235, U-238 and Pu-239.

.. _fig-PE:

.. figure:: _static/images/preequilibrium.pdf
  :width: 50%
  :align: center

  Pre-equilibrium fractions calculated with the :program:`CoH` code. There is only a slight dependence on the target nucleus, and the fit formula (solid line) is used by default in :program:`CGMF` instead.


Excitation Energy, Spin and Parity Distributions
================================================

The total excitation energy (:math:`TXE`) available to the two fragments is constrained by the energy conservation rule

.. math::
   :nowrap:
   :label: TXE

	\begin{eqnarray}
  	TXE &=& Q_f - TKE, \\
      	&=& E_{inc}+B_n+M_n(A_f,Z_f)c^2 - M_n(A_1,Z_1)c^2 - M_n(A_2,Z_2)c^2 - TKE \nonumber
  	\end{eqnarray}

where :math:`TKE` is the total kinetic energy, i.e. the sum of the kinetic energies of fragment 1 and fragment 2, and :math:`M_n` are the nuclear masses for the fissioning nucleus, and the fragments 1 and 2 respectively. Once :math:`TKE` is known, the total excitation energy :math:`TXE` is also known. However, the partitioning of this energy between the two fragments is a more complicated matter, which is discussed at more length in the section below.


Excitation Energy Partitioning
------------------------------

As mentioned above, the total excitation energy (:math:`TXE`) is known as long as the total kinetic energy (:math:`TKE`) and nuclear masses are known. What is not completely known however is the way :math:`TXE` is distributed among the light and the heavy fragments.

Several interesting and competing ideas have been proposed to explain how :math:`TXE` is shared among the two fragments [Schmidt:2010]_ [Talou:2011]_, but no fully compelling proof has been given so far supporting those theories. They all rely on some assumptions regarding the configurations of the fission fragments near the scission point. In the present version of :program:`CGMF`, this excitation energy partitioning is treated as a free parameter, which can be tuned to be best reproduce the average prompt fission neutron multiplicity as a function of the fragment mass, :math:`\overline{\nu}_p(A)`. Indeed, to the first order, the neutron multiplicity reflects the excitation energy of the fragment, while the average neutron energy reflects the temperature of the fragment.

We introduce the ratio of the temperatures between the light and heavy fragments:

.. math::
  :label: RT

  R_T=\frac{T_l}{T_h},

and use the Fermi gas formula to infer the sharing of the excitation energy. This ratio parameter depends on the fragment pair masses :math:`A_l` and :math:`A_h`. At this stage, it is only a convenient way to parameterize the partitioning of :math:`TXE`, and nothing more. Note that this parameter can also be confusing as it uses a ratio of temperatures, while its correct purpose is to share excitation energies. It was introduced at first in the context of the Los Alamos model (LAM) [Madland:1982]_ to compute the average prompt fission neutron spectrum. In its original formulation, the LAM uses a distribution of temperatures to represent the intrinsic excitations in the fragments, and uses the same distribution for both the light and the heavy fragments. In other words, :math:`R_T=1.0`. 

In :program:`CGMF`, :math:`R_T` can be chosen to be mass-dependent to best reproduce :math:`\overline{\nu}_p(A)`. In most cases, it means that :math:`R_T>1.0` as more excitation energy is pumped into the light fragment at the expense of the heavy fragment. This result is in large part due to the deformation energies of the nascent fragments, the heavy fragment being closer to a sphere thanks to shell closures, while the light fragment is largely deformed. This is not true everywhere however, and for very asymmetric fragmentations the inverse becomes true.

We are working on a more physically and mathematically sound proof of this empirical result, in particular in order to expand :program:`CGMF` calculations to other isotopes and energies more reliably.

Figure fig-Ui_ shows an example of a distribution of initial excitation energies in the light and heavy fragments, as well as the total energy, in the case of Cf-252 spontaneous fission.

.. _fig-Ui:

.. figure:: _static/images/Cf252sf_Ui.png
   :width:  70%
   :align:  center

   Typical initial excitation energy distributions in the light and heavy fragments, as well as the total, computed in the case of Cf-252 spontaneous fission.


Spin and Parity Distributions
-----------------------------

The spin of the fragments also follows a conservation rule

.. math::
  :label: spin

  \vec{J_1}+\vec{J_2}+\vec{l}=\vec{J_f}

where :math:`\vec{J_1}` and :math:`\vec{J_2}` are the fission fragment total spins, :math:`\vec{J}` is the total angular momentum of the fissioning nucleus, and :math:`\vec{l}` is the relative orbital angular momentum between the two fragments. In the present version of :program:`CGMF`, :math:`\vec{J_1}` and :math:`\vec{J_2}` follow a Gaussian distribution around a mean value that is chosen to best reproduce some of the observed prompt photon characteristics. The relative orbital angular momentum :math:`l` is left free, so there is no correlation between :math:`\vec{J_1}` and :math:`\vec{J_2}` at this point. This question will be revisited in future versions of the code. Also, negative and positive parities are chosen to be equally probable, so the spin and parity distribution in the fragments reads

.. math::
  :label: JpiDistribution

  \rho(J,\pi) = \frac{1}{2}(2J+1) \exp \left[ -\frac{J(J+1)}{2B^2(Z,A,T)} \right]

where :math:`B` is defined in terms of the fragment temperature as

.. math::

   B^2(Z,A,T)=\alpha\frac{\mathcal{I}_0(A,Z)T}{\hbar^2},

and :math:`\mathcal{I}_0(A,Z)` is the ground-state moment of inertia of the fragment :math:`(A,Z)`. :math:`\alpha` is an adjustable parameter that is used globally to reproduce prompt fission :math:`\gamma` data.

Typical values calculated for the light and heavy fragments are 6-8 :math:`\hbar`, in rather good agreement with values cited in the literature (see [Wilhelmy:1972]_ for instance).



Statistical Hauser-Feshbach Theory
==================================


The Hauser-Feshbach theory [Hauser-Feshbach:1952]_ describes the decay of a compound nucleus in statistical equilibrium through the evaporation of particles and photons until a ground-state or long-lived isomer is reached. This is schematically represented in Fig. 2.6.

.. _fig_diagram:

.. figure:: _static/images/decay_diagram.png
   :width:  70%
   :align:  center

   Schematic drawing explaining the representation of a nucleus in the :program:`CGMF` code, and individual decay paths followed through Monte Carlo simulations.

In this schema, a fragment :math:`(A,Z)` is represented by its ground-state at energy zero, a set of low-lying discrete excited states, and by a set of energy-bins at higher excitation energy where the density of levels becomes too high for individual levels to be separated experimentally. In practice, this picture is not a clear-cut between resolved and unresolved levels. Some levels may have been identified above the continuum threshold region, but it may also be known, from a statistical analysis of the observed levels, that a significant portion of levels has not been observed or that a large fraction of observed levels could not be assigned a specific spin or/and parity. In this case, the matching energy between the discrete and continuum regions is often lowered to well-known levels. 

Fission fragments are neutron-rich, and often relatively far from the valley of :math:`beta`-stability where most experiments have been performed. The known spectroscopy of neutron-rich nuclei is very poor compared to stable nuclei, which means that often very few discrete levels are known. In this case, the matching of the discrete region to the continuum is complicated and very sensitive to the number of specific levels included in the analysis. One also has to rely on systematics of level density parameters to describe the continuum region. Those systematics have been established for stable nuclei and large uncertainties can be expected in the description of nuclei far from stability.

In Fig. 2.6, a couple of decay paths, starting from the same initial excitation energy-bin, are drawn (red arrows) to illustrate the emission of neutrons and photons. In a traditional deterministic Hauser-Feshbach reaction code, the daughter nuclei are all populated at the same time. In a Monte Carlo code such as :program:`CGMF`, only one path is chosen at a given step.

The Hauser-Feshbach theory is statistical in nature and the decay paths are governed by the probabilities for the system to evolve in a particular reaction channel that is open, i.e. physically possible given constraints in energy, spin and parity. We will denote a channel :math:`c` by:

.. math::

	c \equiv (A_i,Z_i,U_i,J_i,\pi_i;A_f,Z_f,U_f,J_f,\pi_f)

In the case of neutron or photon emissions only, we always have :math:`Z_i=Z_f`, and :math:`A_i=A_f` (photon) or :math:`A_f=A_i-1` (neutron).

The probability of decaying through a particular channel :math:`c` is given by the product of the channel transmission coefficients and the density of levels in the final state. For photons, we have:

.. math::

	P(\epsilon_\gamma) dE \propto T_\gamma(\epsilon_\gamma) \rho(Z,A,E-\epsilon_\gamma)dE,

and for neutrons

.. math::

	P(\epsilon_n) dE \propto T_n(\epsilon_n) \rho(Z,A-1,E-\epsilon_n-S_n)dE,

where :math:`\epsilon_\gamma` and :math:`\epsilon_n` are the center-of-mass energies of the emitted photon and neutron, respectively.



Neutron Transmission Coefficients
=================================

Neutron transmission coefficients :math:`T_n^{lj}(\epsilon)` are obtained through optical model calculations. In this model, the Schroedinger equation describing the interaction of incoming waves with a complex mean-field potential is solved, providing the total, shape elastic and reaction cross-sections. It also provides the transmission coefficients that are used in the compound nucleus evaporation calculations.

The transmission coefficients for a channel :math:`c` are obtained from the scattering matrix :math:`S` as

.. math::
  :label: Tn

	T_c=1-\left|\langle S_{cc}\rangle \right|^2.

To calculate the neutron transmission coefficients for fission fragments, it is important to rely on a global optical model potential (OMP) that can provide results for all nuclei. By default, :program:`CGMF` uses the global spherical OMP of Koning and Delaroche [KD03]_.

It is important to note that the calculated spectrum of prompt neutrons does depend on the choice of the optical potential used to compute the neutron transmission coefficients. The OMP of Koning-Delaroche has been established to describe a host of experimental data, e.g., total cross-sections, :math:`S_0` and :math:`S_1` strength functions, etc. However, those data are only available for nuclei near the valley of stability. Some experimental information do indicate that this optical potential may not be very suitable to the fission fragment region, and therefore a relatively large source of uncertainty in the calculation of the neutron spectrum results from this open question.


Gamma-Ray Transmission Coefficients
===================================

The gamma-ray transmission coefficients are obtained using the strength function formalism from the expression: 

.. math::
  :label: Tg

	T^{Xl}(\epsilon_\gamma) = 2\pi f_{Xl}(\epsilon_\gamma)\epsilon_\gamma^{2l+1},

where :math:`\epsilon_\gamma` is the energy of the emitted gamma ray, :math:`Xl` is the multipolarity of the gamma ray, and :math:`f_{Xl}(\epsilon_\gamma)` is the energy-dependent gamma-ray strength function.

For :math:`E1` transitions, the Kopecky-Uhl [Kopecky:1990]_ generalized Lorentzian form for the strength function is used:

.. math::
  :label: E1SF

	f_{E1}(\epsilon_\gamma,T) = K_{E1}\left[ \frac{\epsilon_\gamma \Gamma_{E1}(\epsilon_\gamma)}{\left( \epsilon_\gamma^2-E_{E1}^2\right)^2 + \epsilon^2_\gamma\Gamma_{E1}(\epsilon_\gamma)^2} +\frac{0.7\Gamma_{E1}4\pi^2T^2}{E_{E1}^5} \right] \sigma_{E1}\Gamma_{E1}


where :math:`\sigma_{E1}`, :math:`\Gamma_{E1}`, and :math:`E_{E1}` are the standard giant dipole resonance (GDR) parameters. :math:`\Gamma_{E1}(\epsilon_\gamma)` is an energy-dependent damping width given by

.. math::

	\Gamma_{E1}(\epsilon_\gamma) = \Gamma\frac{\epsilon_\gamma^2+4\pi^2T^2}{E_{E1}^2},

and :math:`T` is the nuclear temperature given by

.. math::

	T=\sqrt{\frac{E^*-\epsilon_\gamma}{a(S_n)}}.

The quantity :math:`S_n` is the neutron separation energy, :math:`E^*` is the excitation energy of the nucleus, and :math:`a` is the level density parameter. The quantity :math:`K_{E1}` is obtained from normalization to experimental data on :math:`2\pi\langle \Gamma_{\gamma_0} \rangle / \langle D_0 \rangle`. 

For :math:`E2` and :math:`M1` transitions, the Brink-Axel [Brink:1955]_ [Axel:1962]_ standard Lorentzian is used instead:

.. math::
  :label: E2M1SF

	f_{Xl}(\epsilon_\gamma)=K_{Xl}\frac{\sigma_{Xl}\epsilon_\gamma\Gamma_{Xl}^2}{(\epsilon_\gamma^2-E_{Xl}^2)^2+\epsilon_\gamma^2\Gamma_{Xl}^2}.


In the current version of :program:`CGMF` (ver. |version|), only :math:`E1, E2`, and :math:`M1` transitions are allowed, and higher multipolarity transitions are neglected.



Level density in the continuum
==============================

In :program:`CGMF`, the Gilbert-Cameron [Gilbert-Cameron:1965]_ model of level densities is used for all fragments. In this model, a constant temperature formula is used to represent the level density at lower excitation energies, while a Fermi gas formula is used at higher excitation energies. Experimental data on the average level spacing at the neutron separation energy can be used to constrain parameters entering the Fermi gas formula, while low-lying discrete levels are used to constrain the constant-temperature parameters. Again, little data is available for nuclei far from stability where systematics have been developed, contribution to uncertainties in the final predicted data.

The constant temperature form is given by

.. math::
  :label: ConstantTemperature

	\rho_{CT}(U)=\frac{1}{T}{\rm exp}\left( \frac{U+\Delta-E_0}{T} \right),

where :math:`T` is the nuclear temperature and :math:`E_0` is a normalization factor. The quantity :math:`U` is the excitation energy :math:`E` minus the pairing energy :math:`\Delta`. At higher excitation energies, the Fermi gas form of the level density is used instead and is given by

.. math::
  :label: FermiGas

	\rho_{FG}(U)=\frac{{\rm exp}\left( 2\sqrt{aU}\right)}{12\sqrt{2}\sigma(U)U(aU)^{1/4}},

where :math:`a` is the level density parameter. The constant temperature form of the level density is matched to cumulative low-lying discrete levels, when they are known. For fission fragments, which are neutron-rich and rather poorly known, this constant-temperature level density is sometimes used down to the ground-state, as shown in the following figure

.. _fig_LD:

.. figure:: _static/images/Rh115_LD.png
   :width:  70%
   :align:  center

In its original formulation, the Gilbert-Cameron formalism uses an energy-independent level density parameter :math:`a`. To better describe the washing-out of shell effects at higher excitation energies, Ignatyuk [Ignatyuk:1979]_ developed a model that uses an energy functional for the level density parameter as

.. math::
  :label: Ignatyuk

	a(U) = \tilde{a} \left( 1+\delta W \frac{1-{\rm exp}(-\gamma U)}{U} \right).

In this formula, :math:`\tilde{a}` is the asymptotic value of the level density parameter at high energy, :math:`\delta W` is the shell correction energy, and :math:`\gamma` is an empirical damping width to account for the washing-out of shell effects at high energy.


Isomeric States
===============

Many low-lying discrete levels that are reported in the ENSDF database have a measurable half-life, ranging from nanoseconds to seconds and even longer. :program:`CGMF` takes this into account when calculating the gamma cascades in the fission products, and samples the exponential decay law according to the reported half-lives. 

An experimental time coincidence window can be set in the ``config.h`` configuration file::

	const double EXPERIMENTAL_TIME_WINDOW = 1e-8;

The time is given in seconds, so in the example above, 1e-8 corresponds to 10 ns. The default value is negative. In this case, all levels are set to decay to the ground-state, ignoring half-lives entirely. Since this value is stored in a configuration file, it is set at compilation time. If the user decides to change this value, he/she would need to recompile the code before using it.

As an example, the calculated intensities for specific gamma lines in Te-134, in the thermal neutron-induced fission of U-235, are shown in the figure below. Time-coincidence windows of 10, 100 and 300 ns were used in three separate calculations. Because of the presence of ~100 ns isomers in Te-134, some of these lines are more or less prominent depending on their half-lives. For example, the :math:`6^+` state at 1.691 MeV has a half-life of 164 ns, decaying to the :math:`4^+` state at 1.576 MeV. A too-short time gate (e.g., 10ns) cannot record this particular gamma line at 115 keV. Similarly, the decay of the :math:`4^+` to :math:`2^+` (297 keV) is also hindered since it depends on the decay of the higher excited :math:`6^+` state.

.. image:: _static/images/134Te.png
   :width: 45%
   :align: left
.. image:: _static/images/isomers.png
   :width: 45%
   :align: right




————————————————————————————————————————————————————————————————————————————————————————————

.. [Randrup:2011] "Brownian Shape Motion on Five-Dimensional Potential-Energy Surfaces: Nuclear Fission-Fragment Mass Distributions," J.Randrup and P.Moller, `Phys. Rev. Lett. 106, 132503 (2011) <http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.132503>`_.

.. [Younes:2014] "A Microscopic Theory of Low Energy Fission: Fragment Properties," W.Younes, D.Gogny, and N.Schunck, in proceedings of the  Fifth International Conference on on ICFN5, "Fission and Properties of Neutron-Rich Nuclei," Sanibel Island, Florida, USA, 4-10 Nov. 2012, Eds. J.H.Hamilton and A.V.Ramaya, World Scientific, p. 605 (2014).

.. [Talou:2011] "Advanced Monte Carlo modeling of prompt fission neutrons for thermal and fast neutron-induced fission reactions on Pu-239," P.Talou, B.Becker, T.Kawano, M.B.Chadwick, and Y.Danon, `Phys. Rev. C 83, 064612 (2011) <http://journals.aps.org/prc/abstract/10.1103/PhysRevC.83.064612>`_.

.. [Wahl:2002] A.C.Wahl, Los Alamos Technical Report LA-13298 (2002).

.. [Madland:2006] "Total prompt energy release in the neutron-induced fission of U-235, U-238, and Pu-239," D. G. Madland, Nucl. Phys. A 772 (2006) 113. 

.. [Viola:1985] V. E. Viola, K. Kwiatkowski, and M. Walker, Phys. Rev. C 31, 1550 (1985).

.. [Hambsch:2007] F.-J. Hambsch (private communication, 2007).

.. [Romano:2010] "Fission fragment mass and energy distributions as a function of incident neutron energy measured in a lead slowing-down spectrometer," C.Romano, Y.Danon, R.Block, J.Thompson, E.Blain, and E.Bond, `Phys. Rev. C 81, 014607 (2010) <http://journals.aps.org/prc/abstract/10.1103/PhysRevC.81.014607>`_.

.. [Wilhelmy:1972] "Angular Momentum of Primary Products Formed in the Spontaneous Fission of Cf-252," J.B.Wilhelmy, E.Cheifetz, R.C.Jared, S.G.Thompson, H.R.Bowman, and J.O.Rasmussen, `Phys. Rev. C 5, 2041 (1972) <http://journals.aps.org/prc/abstract/10.1103/PhysRevC.5.2041>`_.

.. [Schmidt:2010] "Entropy Driven Excitation Energy Sorting in Superfluid Fission Dynamics," K.-H. Schmidt and B. Jurado, `Phys. Rev. Lett. 104, 212501 (2010) <http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.212501>`_.

.. [Hauser-Feshbach:1952] "The inelastic scattering of neutrons," W.Hauser and H.Feshbach, `Phys. Rev. 87, 366-373 (1952) <http://journals.aps.org/pr/abstract/10.1103/PhysRev.87.366>`_. 

.. [Madland:1982] "New calculation of prompt fission neutron spectra and average prompt neutron multiplicities," D.G.Madland and J.R.Nix, Nucl. Sci. Eng. **81**, 213-271 (1982).

.. [KD03] “Local and global nucleon optical models from 1 keV to 200 MeV,” A.J.Koning and J.P.Delaroche, Nucl. Phys. **A713**, 231-310 (2003).

.. [Kopecky:1990] "Test of Gamma-Ray Strength Functions in Nuclear Reaction Model Calculations," J.Kopecky and M.Uhl, `Phys. Rev. C 41, 1941 (1990) <http://journals.aps.org/prc/abstract/10.1103/PhysRevC.41.1941>`_.

.. [Brink:1955] "Some aspects of the interaction of fields with matter," D.M.Brink, D. Ph. Thesis, Oxford (1955); "Individual Particle and Collective Aspects of the Nuclear Photoeffect," Nucl. Phys. **4**, 215 (1957).

.. [Axel:1962] "Electric Dipole Ground State Transition Width Strength Function and 7-MeV Photon Interactions," P.Axel, `Phys. Rev. 126, 671 (1962) <http://journals.aps.org/pr/abstract/10.1103/PhysRev.126.671>`_.

.. [Gilbert-Cameron:1965] "A composite nuclear level density formula with shell corrections," A.Gilbert and A.G.W.Cameron, Can. J. Phys. **43**, 1446 (1965).

.. [Ignatyuk:1979] A.V.Ignatyuk, K.K.Istekov, and G.N.Smirenkin, Sov. J. Nucl. Phys. **29**, 450 (1979).


