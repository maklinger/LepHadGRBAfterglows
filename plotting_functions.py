import numpy as np
from astropy.cosmology import Planck18

# for absorption ranges with fading
def plot_logfaderange(ax, E1, E2, facE1, facE2, color, N=5):
    # main range alpha=1
    ax.axvspan(E1, E2, color=color, alpha=1)
    # N bins to fade out
    bounds1 = np.logspace(np.log10(E1/facE1), np.log10(E1), N)
    bounds2 = np.logspace(np.log10(E2), np.log10(E2*facE2), N)
    for i in range(N-1):
        ax.axvspan(bounds1[i], bounds1[i+1], color=color, alpha=(i+1)/N)
        ax.axvspan(bounds2[i], bounds2[i+1], color=color, alpha=1-(i+1)/N)
    return ax


def dens2flux_av(EdN_dlgE, t_obs, z, Gamma):
    """
    converts density into flux [erg/s/cm^2] for Gamma>>1 with average angle cos(theta) = beta
    args:
        EdN_dlgE: comov. number spectrum of free streaming particles [erg]
        t_obs   : observed time [s]
        z       : redshift
        Gamma   : Lorentz factor of shell
    returns:
        observed energy flux [erg/(cm^2s)]
    """
    c = 2.99792458e10 # cm/s
    r = 4*Gamma**2 * t_obs * c
    dL = Planck18.luminosity_distance(z).to("cm").value
    return EdN_dlgE * r**3/t_obs / dL**2


def Eobs_av(E_co, z=0, Gamma=1):
    """
    converts comoving to observed energy with average angle cos(theta) = beta
    args:
        E_co   : comoving energy [any energy unit]
        z      : redshift
        Gamma  : Lorentz factor of blob
    returns:
        observed energy [same unit as E_co]
    """
    return E_co * Gamma / (1+z)


def plot_spectrum(
    ax, z, Gamma, t_obs, Eg_eV, En_eV,
    E2dN_dE_tot_gamma, E2dN_dE_tot_nu, components,
    cgrey="#444444", Emins_grey=[1, 1e13], Emaxs_grey=[3e2, 1e22], facE_grey=5, Ngrey=30,
    labels_grey=["dust\nphotoel.", r"$\gamma \gamma$ (EBL,CMB)"], textcolgrey="lightgrey",
    Emins_obs=[1e3, 3e11], Emaxs_obs=[1e5, 1e13],
    dlgobs=-0.7, color_obs="lightgrey", labels_obs=["X-ray", "VHE"],
    alpha_dom=0.4, zorder_dom=10, lw_main=3, col_tot="k", col_nu="teal",
    ylim=[1e-18, 1e-5], xlim=[1e-4, 1e18], tobs="1000", Ekiniso=None
):

    # grey band for absorption
    ygrey = 0.8 * ylim[1]
    for i in range(len(Emins_grey)):
        ax = plot_logfaderange(
            ax, Emins_grey[i], Emaxs_grey[i], facE_grey, facE_grey, cgrey, Ngrey)
        # label for grey band
        ax.text(np.sqrt(Emins_grey[i]*Emaxs_grey[i]), ygrey, labels_grey[i],
                va="top", ha="center", ma="center", color=textcolgrey)

    # observational ranges
    ymax_obs = ylim[1]
    ymin_obs = 10**dlgobs * ymax_obs
    for i in range(len(Emins_obs)):
        ax.fill_between([Emins_obs[i], Emaxs_obs[i]], ymin_obs*np.ones(2),
                        ymax_obs*np.ones(2), alpha=1, color=color_obs)
        ax.text(np.sqrt(Emins_obs[i]*Emaxs_obs[i]), np.sqrt(ymin_obs*ymax_obs),
                labels_obs[i], va="center", ha="center", ma="center", color="k")

    Egobs = Eobs_av(Eg_eV, z, Gamma)

    # total photons
    ax.loglog(
        Egobs, dens2flux_av(E2dN_dE_tot_gamma, t_obs, z, Gamma), lw=lw_main,
        zorder=zorder_dom, color=col_tot)

    # total neutrinos
    ax.loglog(
        Eobs_av(En_eV, z, Gamma),
        dens2flux_av(E2dN_dE_tot_nu, t_obs, z, Gamma),
        lw=lw_main, zorder=zorder_dom, color=col_nu)

    # components

    for E2dN_dE, col, dom, zup, ls, lbl, lw in components:
        if dom:
            ax.fill_between(
                Egobs, np.zeros_like(Egobs),
                dens2flux_av(E2dN_dE, t_obs, z, Gamma),
                color=col, alpha=alpha_dom, zorder=zorder_dom-1+zup, ls=ls,
                lw=lw)
        ax.loglog(
            Egobs, dens2flux_av(E2dN_dE, t_obs, z, Gamma),
            color=col, ls=ls, zorder=zorder_dom+zup, label=lbl, lw=lw)

    title_string = r"$t=" + f"{tobs}"
    title_string += r"$ s, $z=%g" % (z)
    if Ekiniso is not None:
        title_string += r",\; E_\mathrm{kin,iso}=%.0e \; \mathrm{erg}" %(Ekiniso)
    title_string += "$"
    ax.set_title(title_string, loc="left")
    ax.set_ylim(*ylim)
    ax.set_yticks(10**np.arange(np.log10(ax.get_ylim()
                  [0]), np.log10(ax.get_ylim()[1])+0.5))
    ax.set_xlim(*xlim)
    ax.set_xticks(10**np.arange(-3, np.log10(ax.get_xlim()[1])+0.5, 3))
    ax.set_xticks(10**np.arange(np.log10(ax.get_xlim()
                  [0]), np.log10(ax.get_xlim()[1])+0.5, 1), minor=True)
    ax.set_xticklabels(["" for i in np.arange(
        np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1])+0.5, 1)], minor=True)
    ax.set_aspect("equal")
    ax.grid(which="both", alpha=0.3, zorder=-10)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$E$ [eV]")
    ax.set_ylabel(r"$EF_E$ [erg/(cm²s)]")
    return ax



def plot_coolingtimes(
        ax, Ep_eV, Ee_eV, Eg_eV, ind, Eemin_eV, Epmin_eV,
        tpsys, tpads, tpacs, tpics, tpbhs, tppgs, tppps,
        tesys, teics, teacs, tgescs, tgpairs, tgssas, tpidec, tmudec,
        cel="tab:blue", cproton="tab:purple", cgamma="tab:green", cpion="tab:cyan",
        cmuon="tab:grey", cadi="k", alpha_inj=0.2,
        ylim=[1e-1, 1e15], xlim=[1e-4, 1e21], leptonic=False,
        pions=True, muons=True):

    # injection ranges
    tels = 1/(1/tesys[ind] + 1/tpads[ind][0] + 1/teics[ind])
    Eemax_eV = Ee_eV[np.argmin((tels-teacs[ind])**2)]
    tps = 1/(1/tpsys[ind] + 1/tpads[ind][0] + 1/tpics[ind] +
             1/tpbhs[ind] + 1/tppgs[ind] + 1/tppps[ind])
    Epmax_eV = Ep_eV[np.argmin((tps-tpacs[ind])**2)]
    if not leptonic:
        ax.axvspan(Epmin_eV, Epmax_eV, color=cproton, alpha=alpha_inj)
        ax.axvline(Epmax_eV, color=cproton)
    ax.axvspan(Eemin_eV, Eemax_eV, color=cel, alpha=alpha_inj)
    ax.axvline(Eemax_eV, color=cel)
    ax.axvline(Eemin_eV, color=cel)
    if not leptonic:
        ax.axvline(Epmin_eV, color=cproton, ls="--")

    if not leptonic:
        massp = 0.938  # proton mass / GeV
        masspipm = 0.1396
        massmu = 0.10566
        if muons:
            ax.loglog(Ep_eV, (massmu/massp)**4 *
                    tpsys[ind], label="mu syn", c=cmuon)
            ax.loglog(Ep_eV, tpidec, label="pi dec", c=cpion, ls=":")
        if pions:
            ax.loglog(Ep_eV, (masspipm/massp)**4 *
                    tpsys[ind], label="pi syn", c=cpion)
            ax.loglog(Ep_eV, tmudec, label="mu dec", c=cmuon, ls=":")

    # esc
    ax.axhline(tgescs[ind][0], label="g esc", c=cgamma)

    # adi
    if leptonic:
        cadi = cel
    ax.axhline(tpads[ind][0], label="e adi", c=cadi, ls="-")
    # if not leptonic:
    #     ax.axhline(tpads[ind][0], label="p adi", c=cproton, ls="--")

    # acc
    ax.loglog(Ee_eV, teacs[ind], label="e acc", c=cel, ls="-")
    if not leptonic:
        ax.loglog(Ep_eV, tpacs[ind], label="p acc", c=cproton, ls="--")

    ax.loglog(Ee_eV, tesys[ind], label="e syn", c=cel, ls="-")
    ax.loglog(Ee_eV, teics[ind], label="e ic", c=cel, ls="--")

    if not leptonic:
        ax.loglog(Ep_eV, tpsys[ind], label="p syn", c=cproton)
        ax.loglog(Ep_eV, tpics[ind], label="p ic", c=cproton, ls="--")
        ax.loglog(Ep_eV, tpbhs[ind], label="p bh", c=cproton, ls=":")
        ax.loglog(Ep_eV, tppgs[ind], label="p pg", c=cproton, ls="-.")
        ax.loglog(Ep_eV, tppps[ind], label="p pp", c=cproton, ls="-")

    ax.loglog(Eg_eV, tgpairs[ind], label="g pair", c=cgamma, ls="--")
    ax.loglog(Eg_eV, tgssas[ind], label="p ssa", c=cgamma, ls="-.")

    ax.set_ylim(*ylim)
    ax.set_yticks(10**np.arange(np.log10(ax.get_ylim()
                  [0]), np.log10(ax.get_ylim()[1])+0.5))
    ax.set_xlim(*xlim)
    ax.set_xticks(10**np.arange(-3, np.log10(ax.get_xlim()[1])+0.5, 3))
    ax.set_xticks(10**np.arange(np.log10(ax.get_xlim()
                  [0]), np.log10(ax.get_xlim()[1])+0.5, 1), minor=True)
    ax.set_xticklabels(["" for i in np.arange(
        np.log10(ax.get_xlim()[0]), np.log10(ax.get_xlim()[1])+0.5, 1)], minor=True)
    ax.set_aspect("equal")
    ax.grid(which="both", alpha=0.3, zorder=-10)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$E'$ [eV]")
    ax.set_ylabel(r"$\tau'$ [s]")
    return ax
