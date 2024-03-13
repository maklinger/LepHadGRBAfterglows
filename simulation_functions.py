import am3lib
import numpy as np
import physicsConsts as c

# this defines the standard set of switches
def get_AM3(exp=1, hadronic=1, qsyn=0):
    am3 = am3lib.AM3()
    am3.update_energy_grid(1e-6, 1e3, 1e22)
    # set switches
    am3.set_estimate_max_energies(0)

    am3.set_process_parse_sed(1)
    am3.set_process_hadronic(hadronic)
    am3.set_process_merge_positrons_into_electrons(0)
    am3.set_process_escape(1)
    # expansion related
    am3.set_process_adiabatic_cooling(1)
    am3.set_process_expansion(exp)

    # synchrotron related
    am3.set_process_electron_syn(1)
    am3.set_process_electron_syn_emission(1)
    am3.set_process_electron_syn_cooling(1)
    am3.set_process_quantum_syn(qsyn)
    am3.set_process_ssa(1)
    am3.set_process_proton_syn(1)
    am3.set_process_proton_syn_emission(1)
    am3.set_process_proton_syn_cooling(1)
    am3.set_process_pion_syn(1)
    am3.set_process_pion_syn_emission(1)
    am3.set_process_pion_syn_cooling(1)
    am3.set_process_muon_syn(1)
    am3.set_process_muon_syn_emission(1)
    am3.set_process_muon_syn_cooling(1)

    # inverse Compton related
    am3.set_process_electron_compton(1)
    am3.set_process_electron_compton_emission(1)
    am3.set_process_electron_compton_cooling(1)
    am3.set_process_compton_photon_energy_loss(1)
    am3.set_process_proton_compton(1)
    am3.set_process_proton_compton_emission(1)
    am3.set_process_proton_compton_cooling(1)
    am3.set_process_pion_compton(1)
    am3.set_process_pion_compton_emission(1)
    am3.set_process_pion_compton_cooling(1)
    am3.set_process_muon_compton(1)
    am3.set_process_muon_compton_emission(1)
    am3.set_process_muon_compton_cooling(1)

    # secondary decay
    am3.set_process_pion_decay(1)
    am3.set_process_muon_decay(1)

    # pair production (gamma+gamma -> e- + e+)
    am3.set_process_annihilation(1)
    am3.set_process_annihilation_cooling(1)
    am3.set_process_annihilation_pair_emission(1)

    # Bethe-Heitler
    am3.set_process_bethe_heitler(1)
    am3.set_process_bethe_heitler_emission(1)
    am3.set_process_bethe_heitler_cooling(1)

    # p-gamma
    am3.set_process_photopion(1)
    am3.set_process_photopion_emission(1)
    am3.set_process_photopion_cooling(1)
    am3.set_process_photopion_photon_loss(1)

    # proton proton
    am3.set_process_pp(1)
    am3.set_process_pp_emission(1)
    am3.set_process_pp_emission_pi0_into_cascade(0)
    am3.set_process_pp_cooling(1)


    # optimisations

    am3.set_optimize_compton_emission_emin(1e-3)
    am3.set_optimize_compton_emission_grid(1)
    am3.set_optimize_compton_target_emax(1e6)
    am3.set_optimize_compton_target_grid(1)

    am3.set_optimize_annihilation_pair_emission(1)

    am3.set_optimize_bethe_heitler_outgoing_pairs_grid(1)
    am3.set_optimize_bethe_heitler_incoming_protons_min(1e12)
    am3.set_optimize_bethe_heitler_target_photon_max(1e6)

    am3.set_optimize_photopion_target_photon_grid(1)
    am3.set_optimize_photopion_target_photon_max(1e6)

    am3.set_pp_model_charged_pions(0)
    am3.set_pp_model_neutral_pions(0)

    am3.set_pp_transition_energy_to_delta_approx(1e11)

    # am3.set_solver_simpleqintegration(solver_simpleqintegration)
    am3.set_solver_threshold_cool_dom(1e-1)
    # am3.set_solver_threshold_esc_dom(solver_threshold_esc_dom)
    # am3.set_solver_threshold_matrix(solver_threshold_matrix)

    return am3


def powerlaw_array(Es, Emin, Emax, p, norm):
    powerlaw = (Es/Emin)**(-p) * np.exp(-(Es/Emax)) * np.greater_equal(Es, Emin)
    integral = np.log(Es[1]/Es[0]) * np.sum(Es**2 * powerlaw)
    return norm/integral * powerlaw

def Gamma_tobs(t_obs, G0, t0):
    return G0 * (t_obs/t0)**(-3/8)

def blastwave(t_eval, Gamma, n, t_obs, epsilon_B):
    Gamma_eval = Gamma * (t_eval/t_obs)**(-3/8)
    t_dyn = Gamma_eval * t_eval
    n_target = Gamma_eval * n
    p_ram_eval = n * Gamma_eval**2 * c.mpc2
    B = np.sqrt(8*np.pi*epsilon_B * p_ram_eval)
    return t_dyn, n_target, p_ram_eval, B