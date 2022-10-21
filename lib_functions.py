import math
import numpy as np
import matplotlib.markers as mark
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Helvetica']})
rc('text', usetex=True)

budget = 1
eta_h = 1
eta_l = 0.50
theta = 0.50
delta = 1
beta = 2
alpha = 1
gamma_l = 0.50
prob = np.linspace(.01, .99, 100)
iter = 5000


# donor_alpha, market_beta, npo_delta, market_theta, low_eta, high_eta, prob_shock, low_gamma, donor_budget

def higher_shock(p, low_gamma):
    try:
        gamma_h = (1 - p * low_gamma) / (1 - p)
    except ZeroDivisionError:
        gamma_h = 1
    return gamma_h


def eta_bar(npo_theta, low_eta, high_eta=eta_h):
    return npo_theta * low_eta + (1 - npo_theta) * high_eta


def donor_utility_tf(npo_theta, low_eta, donor_alpha, donor_budget=budget):
    if donor_alpha <= 2 * eta_bar(npo_theta, low_eta):
        u_d = 2 * donor_budget * eta_bar(npo_theta, low_eta)
        y = 2 * donor_budget * eta_bar(npo_theta, low_eta)
    else:
        u_d = donor_alpha * donor_budget
        y = 0
    return u_d, y


def npo_expected_y_tf(low_eta, donor_alpha, npo_theta, high_eta=eta_h, donor_budget=budget):
    if donor_alpha <= 2 * eta_bar(npo_theta, low_eta):
        y_l = 2 * donor_budget * low_eta
        y_h = 2 * donor_budget * high_eta
    else:
        y_l = 0
        y_h = 0
    return y_l, y_h


def u_n_tf(low_eta, donor_alpha, npo_theta, high_eta=eta_h, donor_budget=budget):
    if donor_alpha <= 2 * eta_bar(npo_theta, low_eta):
        u_n_l = donor_budget * low_eta
        u_n_h = donor_budget * high_eta
    else:
        u_n_l = 0
        u_n_h = 0
    return u_n_l, u_n_h


def ratio_to_region(ratio, p, low_gamma):
    high_gamma = higher_shock(p, low_gamma)
    m_b = (2 * high_gamma / (2 * high_gamma - 1)) ** 2
    m_c = (1 / p) * (2 * high_gamma / (2 * high_gamma - low_gamma)) ** 2
    if (ratio >= 1) and (ratio < m_b):
        region = 1
    elif (ratio >= m_b) and (ratio < min(m_c, 4)):
        region = 2
    elif (ratio >= m_c) and (ratio < 4):
        region = 3
    elif (ratio >= 4) and (ratio < m_c):
        region = 4
    elif (ratio >= max(m_c, 4)) and (ratio < 4 / p):
        region = 5
    else:
        region = 6
    return region


def yt_thresholds(market_beta, npo_delta, eta, prob_shock, low_gamma, donor_budget=budget):
    high_gamma = higher_shock(prob_shock, low_gamma)

    y_t_a = donor_budget * eta * low_gamma * (
            1 + math.sqrt(1 + (2 * npo_delta * market_beta * prob_shock / eta))) / market_beta
    y_t_b = donor_budget * eta * high_gamma * \
            (1 + math.sqrt(
                1 + (2 * (low_gamma / high_gamma) * npo_delta * market_beta * prob_shock / eta))) / market_beta
    try:
        y_t_c = 2 * donor_budget * eta * high_gamma * \
                (1 + math.sqrt(
                    1 - ((1 - (low_gamma / high_gamma)) * npo_delta * market_beta * prob_shock / eta))) / market_beta
    except ValueError:
        y_t_c = 0

    try:
        y_t_c_e = 2 * donor_budget * eta * high_gamma * \
                  (1 - math.sqrt(
                      1 - ((1 - (
                              low_gamma / high_gamma)) * npo_delta * market_beta * prob_shock / eta))) / market_beta
    except ValueError:
        y_t_c_e = 0

    y_t_d = npo_delta * donor_budget / ((math.sqrt(npo_delta * market_beta / eta)) - 1)
    return y_t_a, y_t_b, y_t_c, y_t_d, y_t_c_e


def delta_beta_thresholds(prob_shock, low_gamma):
    high_gamma = higher_shock(prob_shock, low_gamma)

    m_a = ((2 * high_gamma) / (3 * high_gamma - prob_shock * high_gamma - 1)) ** 2
    m_b = ((2 * high_gamma) / (2 * high_gamma - 1)) ** 2
    m_c = (4 * prob_shock * high_gamma ** 2) / ((high_gamma * prob_shock + high_gamma - 1) ** 2)
    return m_a, m_b, m_c


# Given a target y_T, the maximum npo's utility and corresponding y delivered.
def find_expected_y_l(yt, p, low_gamma, market_beta, eta, npo_delta=delta, donor_budget=budget):
    high_gamma = higher_shock(p, low_gamma)
    npo_regional_utility = []
    y_region_utility = []
    ratio = npo_delta * market_beta / eta
    y_t_a, y_t_b, y_t_c, y_t_d, y_t_c_e = yt_thresholds(market_beta, npo_delta, eta, p, low_gamma, donor_budget=budget)

    if yt <= 2 * donor_budget * eta * low_gamma / market_beta:
        npo_regional_utility.append(donor_budget * eta / market_beta)
        y_region_utility.append(2 * donor_budget * eta / market_beta)
        if donor_budget * eta / market_beta < 0:
            print("region1")

    if 2 * donor_budget * eta * low_gamma / market_beta < yt < 4 * donor_budget * eta * low_gamma / market_beta:
        npo_regional_utility.append((yt / low_gamma) * (1 - (yt * market_beta / (4 * donor_budget * eta * low_gamma))))
        y_region_utility.append(yt / low_gamma)
        if (yt / low_gamma) * (1 - (yt * market_beta / (4 * donor_budget * eta * low_gamma))) < 0:
            print("region2")

    if (2 * donor_budget * low_gamma / market_beta) * math.sqrt(
            npo_delta * market_beta * eta) < yt <= y_t_d * low_gamma and \
            4 > ratio >= 1:
        npo_regional_utility.append(
            npo_delta * donor_budget + (yt / low_gamma) * (1 - math.sqrt(npo_delta * market_beta / eta)))
        y_region_utility.append(yt / low_gamma)
        if npo_delta * donor_budget + (yt / low_gamma) * (1 - math.sqrt(npo_delta * market_beta / eta)) < 0:
            print("region3")

    if y_t_a < yt <= y_t_b and (4 / p) > ratio > 1 and \
            ratio < (1 / p) * (1 + (npo_delta * donor_budget * p * low_gamma / yt)) ** 2:
        npo_regional_utility.append(
            (donor_budget / market_beta) * eta * (1 + (npo_delta * p * low_gamma * donor_budget / yt)) ** 2
            - npo_delta * donor_budget * p)
        y_region_utility.append(
            2 * donor_budget * eta * (1 + (npo_delta * p * low_gamma * donor_budget / yt)) / market_beta)
        if ((donor_budget / market_beta) * eta * (1 + (npo_delta * p * low_gamma * donor_budget / yt)) ** 2
            - npo_delta * donor_budget * p) < 0:
            print("region4")

    if max(y_t_c_e, y_t_b,
           (2 * donor_budget * high_gamma / market_beta) * math.sqrt(eta * market_beta * p * npo_delta)) < yt <= y_t_c \
            and ratio > 1:
        npo_regional_utility.append((yt / high_gamma) * (1 - (yt * market_beta / (4 * donor_budget * eta * high_gamma)))
                                    - npo_delta * donor_budget * p * (1 - (low_gamma / high_gamma)))
        y_region_utility.append(yt / high_gamma)
        if (yt / high_gamma) * (1 - (yt * market_beta / (4 * donor_budget * eta * high_gamma))) \
                - npo_delta * donor_budget * p * (1 - (low_gamma / high_gamma)) < 0:
            print("region5")

    if 2 * donor_budget * high_gamma * math.sqrt(npo_delta * market_beta * eta) < yt <= y_t_d and \
            (1 + (npo_delta * donor_budget * p * low_gamma / yt)) ** 2 < ratio <= \
            delta_beta_thresholds(p, low_gamma)[1]:
        npo_regional_utility.append((npo_delta * donor_budget / high_gamma) +
                                    (yt / high_gamma) * (1 - math.sqrt(npo_delta * market_beta / eta)))
        y_region_utility.append(yt / high_gamma)
        if (npo_delta * donor_budget / high_gamma) + \
                (yt / high_gamma) * (1 - math.sqrt(npo_delta * market_beta / eta)) < 0:
            print("region6")

    npo_regional_utility.append(0)
    y_region_utility.append(0)

    u_n = max(npo_regional_utility)
    index = npo_regional_utility.index(max(npo_regional_utility))
    y = y_region_utility[index]

    return u_n, y


def donation_pfr(yt, y, p, low_gamma, donor_budget=budget):
    high_gamma = higher_shock(p, low_gamma)
    if y * low_gamma >= yt:
        d = donor_budget
    elif y * low_gamma < yt <= y * high_gamma:
        d = (p * y * low_gamma * donor_budget / yt) + (1 - p) * donor_budget
    else:
        d = y * donor_budget / yt
    return d


def donor_utility_pfr_yt(yt, p, low_gamma, low_eta, donor_alpha, npo_theta, market_beta, npo_delta,
                         donor_budget=budget, high_eta=eta_h):
    y_l = find_expected_y_l(yt, p, low_gamma, market_beta, low_eta, npo_delta, donor_budget=budget)[1]
    u_nl = find_expected_y_l(yt, p, low_gamma, market_beta, low_eta, npo_delta, donor_budget=budget)[0]
    u_nh = find_expected_y_l(yt, p, low_gamma, market_beta, high_eta, npo_delta, donor_budget=budget)[0]
    y_h = find_expected_y_l(yt, p, low_gamma, market_beta, high_eta, npo_delta, donor_budget=budget)[1]
    d_l = donation_pfr(yt, y_l, p, low_gamma)
    d_h = donation_pfr(yt, y_h, p, low_gamma)
    u_d = npo_theta * (y_l + donor_alpha * (donor_budget - d_l)) \
          + (1 - npo_theta) * (y_h + donor_alpha * (donor_budget - d_h))
    expected_y = (npo_theta * y_l) + (1 - npo_theta) * y_h
    return u_d, expected_y, u_nl, u_nh, y_l, y_h


def max_yt_possible(market_beta, npo_delta, p, low_gamma, high_eta=eta_h, donor_budget=budget):
    ratio_h = npo_delta * market_beta / high_eta
    a = max(yt_thresholds(market_beta, npo_delta, high_eta, p, low_gamma, donor_budget=budget))
    b = p * low_gamma * npo_delta * donor_budget / (math.sqrt(ratio_h * p) - 1)
    c = 4 * donor_budget * high_eta * low_gamma / market_beta
    return max([a, b, c])


def plots_with_gamma_l(gamma=np.linspace(.01, .05, 3), prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                       donor_budget=budget,
                       npo_delta=delta, market_beta=beta, donor_alpha=alpha, npo_theta=theta, num=iter):
    donor_pfr_optimal_dict = {}
    donor_pfr_yt_dict = {}
    y_pfr_optimal_dict = {}
    y_tf_optimal_dict = {}
    u_d_tf_optimal_dict = {}
    npo_pfr_optimal_low = {}
    npo_pfr_optimal_high = {}
    range_dic_l = {}
    range_dic_h = {}
    ratio_h = npo_delta * market_beta / high_eta
    ratio_l = npo_delta * market_beta / low_eta
    y_l = {}
    y_h = {}

    for gammaL in gamma:
        u_d_pfr_with_p_dict = {}
        yt_with_p_dict = {}
        u_nl_with_p_dict = {}
        u_nh_with_p_dict = {}
        expected_y_with_p_dict = {}
        range_dict_l_with_p = {}
        range_dict_h_with_p = {}
        y_l_pfr_with_p = {}
        y_h_pfr_with_p = {}
        for p in prob_vector:

            range_dict_l_with_p[p] = ratio_to_region(ratio_l, p, low_gamma=gammaL)
            range_dict_h_with_p[p] = ratio_to_region(ratio_h, p, low_gamma=gammaL)

            bar_yt = 2 * max_yt_possible(market_beta, npo_delta, p, gammaL, high_eta, donor_budget)
            y = np.linspace(0, bar_yt, num)
            u_d_pfr_with_yt_dict = {}

            for j in y:
                u_d_pfr_with_yt_dict[j] = donor_utility_pfr_yt(j, p, gammaL, low_eta, donor_alpha, npo_theta,
                                                               market_beta, npo_delta, donor_budget,
                                                               high_eta)

            yt_optimal = max(u_d_pfr_with_yt_dict, key=lambda k: u_d_pfr_with_yt_dict[k][0])

            expected_y_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][1]
            u_d_pfr_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][0]  # for a given p, this is optimal u_d_pfr
            u_nl_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][2]
            u_nh_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][3]
            y_l_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][4]
            y_h_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][5]
            yt_with_p_dict[p] = yt_optimal

        donor_pfr_optimal_dict[gammaL] = [u_d_pfr_with_p_dict[k] for k in prob_vector]
        y_pfr_optimal_dict[gammaL] = [expected_y_with_p_dict[k] for k in prob_vector]
        donor_pfr_yt_dict[gammaL] = [yt_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_low[gammaL] = [u_nl_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_high[gammaL] = [u_nh_with_p_dict[k] for k in prob_vector]
        y_l[gammaL] = [y_l_pfr_with_p[k] for k in prob_vector]
        y_h[gammaL] = [y_h_pfr_with_p[k] for k in prob_vector]
        u_d_tf_optimal_dict[gammaL] = [donor_utility_tf(npo_theta, low_eta, donor_alpha)[0]] * len(prob_vector)
        y_tf_optimal_dict[gammaL] = [donor_utility_tf(npo_theta, low_eta, donor_alpha)[1]] * len(prob_vector)
        range_dic_l[gammaL] = [range_dict_l_with_p[k] for k in prob_vector]
        range_dic_h[gammaL] = [range_dict_h_with_p[k] for k in prob_vector]
    return donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
           npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h


def plots_with_theta(theta_vector, prob_vector, low_eta, high_eta,
                     donor_budget,
                     npo_delta, market_beta, donor_alpha, low_gamma, num=iter):
    donor_pfr_optimal_dict = {}
    donor_pfr_yt_dict = {}
    u_d_tf_optimal_dict = {}
    y_pfr_optimal_dict = {}
    y_tf_optimal_dict = {}
    npo_pfr_optimal_low = {}
    npo_pfr_optimal_high = {}
    range_dic_l = {}
    range_dic_h = {}
    ratio_h = npo_delta * market_beta / high_eta
    ratio_l = npo_delta * market_beta / low_eta
    y_l = {}
    y_h = {}

    for theta_value in theta_vector:
        u_d_pfr_with_p_dict = {}
        yt_with_p_dict = {}
        u_nl_with_p_dict = {}
        u_nh_with_p_dict = {}
        expected_y_with_p_dict = {}
        range_dict_l_with_p = {}
        range_dict_h_with_p = {}
        y_l_pfr_with_p = {}
        y_h_pfr_with_p = {}
        for p in prob_vector:

            range_dict_l_with_p[p] = ratio_to_region(ratio_l, p, low_gamma)
            range_dict_h_with_p[p] = ratio_to_region(ratio_h, p, low_gamma)
            bar_yt = 2 * max_yt_possible(market_beta, npo_delta, p, low_gamma, high_eta, donor_budget)
            y = np.linspace(0, bar_yt, num)
            u_d_pfr_with_yt_dict = {}

            for j in y:
                u_d_pfr_with_yt_dict[j] = donor_utility_pfr_yt(j, p, low_gamma, low_eta, donor_alpha, theta_value,
                                                               market_beta, npo_delta,
                                                               donor_budget, high_eta)

            yt_optimal = max(u_d_pfr_with_yt_dict, key=lambda k: u_d_pfr_with_yt_dict[k][0])

            expected_y_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][1]
            u_d_pfr_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][0]  # for a given p, this is optimal u_d_pfr
            u_nl_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][2]
            u_nh_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][3]
            y_l_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][4]
            y_h_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][5]
            yt_with_p_dict[p] = yt_optimal

        donor_pfr_optimal_dict[theta_value] = [u_d_pfr_with_p_dict[k] for k in prob_vector]
        donor_pfr_yt_dict[theta_value] = [yt_with_p_dict[k] for k in prob_vector]
        y_pfr_optimal_dict[theta_value] = [expected_y_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_low[theta_value] = [u_nl_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_high[theta_value] = [u_nh_with_p_dict[k] for k in prob_vector]
        y_l[theta_value] = [y_l_pfr_with_p[k] for k in prob_vector]
        y_h[theta_value] = [y_h_pfr_with_p[k] for k in prob_vector]
        u_d_tf_optimal_dict[theta_value] = [donor_utility_tf(theta_value, low_eta, donor_alpha,
                                                             donor_budget=budget)[0]] * len(prob_vector)
        y_tf_optimal_dict[theta_value] = [donor_utility_tf(theta_value, low_eta, donor_alpha,
                                                           donor_budget=budget)[1]] * len(prob_vector)
        range_dic_l[theta_value] = [range_dict_l_with_p[k] for k in prob_vector]
        range_dic_h[theta_value] = [range_dict_h_with_p[k] for k in prob_vector]

    return donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
           npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h


def plots_with_alpha(alpha_vector, prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                     donor_budget=budget,
                     npo_delta=delta, market_beta=beta, npo_theta=theta, low_gamma=gamma_l, num=iter):
    donor_pfr_optimal_dict = {}
    donor_pfr_yt_dict = {}
    u_d_tf_optimal_dict = {}
    y_pfr_optimal_dict = {}
    npo_pfr_optimal_low = {}
    npo_pfr_optimal_high = {}
    y_tf_optimal_dict = {}
    range_dic_l = {}
    range_dic_h = {}
    ratio_h = npo_delta * market_beta / high_eta
    ratio_l = npo_delta * market_beta / low_eta
    y_l = {}
    y_h = {}

    for alpha_value in alpha_vector:
        u_d_pfr_with_p_dict = {}
        yt_with_p_dict = {}
        u_nl_with_p_dict = {}
        u_nh_with_p_dict = {}
        expected_y_with_p_dict = {}
        range_dict_l_with_p = {}
        range_dict_h_with_p = {}
        y_l_pfr_with_p = {}
        y_h_pfr_with_p = {}
        for p in prob_vector:

            range_dict_l_with_p[p] = ratio_to_region(ratio_l, p, low_gamma)
            range_dict_h_with_p[p] = ratio_to_region(ratio_h, p, low_gamma)
            bar_yt = 2 * max_yt_possible(market_beta, npo_delta, p, low_gamma, high_eta, donor_budget)
            y = np.linspace(0, bar_yt, num)

            u_d_pfr_with_yt_dict = {}
            for j in y:
                u_d_pfr_with_yt_dict[j] = donor_utility_pfr_yt(j, p, low_gamma, low_eta, alpha_value, npo_theta,
                                                               market_beta, npo_delta,
                                                               donor_budget, high_eta)

            yt_optimal = max(u_d_pfr_with_yt_dict, key=lambda k: u_d_pfr_with_yt_dict[k][0])

            expected_y_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][1]
            u_d_pfr_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][0]  # for a given p, this is optimal u_d_pfr
            u_nl_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][2]
            u_nh_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][3]
            y_l_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][4]
            y_h_pfr_with_p[p] = u_d_pfr_with_yt_dict[yt_optimal][5]
            yt_with_p_dict[p] = yt_optimal

        donor_pfr_optimal_dict[alpha_value] = [u_d_pfr_with_p_dict[k] for k in prob_vector]
        donor_pfr_yt_dict[alpha_value] = [yt_with_p_dict[k] for k in prob_vector]
        y_pfr_optimal_dict[alpha_value] = [expected_y_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_low[alpha_value] = [u_nl_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_high[alpha_value] = [u_nh_with_p_dict[k] for k in prob_vector]
        y_l[alpha_value] = [y_l_pfr_with_p[k] for k in prob_vector]
        y_h[alpha_value] = [y_h_pfr_with_p[k] for k in prob_vector]
        u_d_tf_optimal_dict[alpha_value] = [donor_utility_tf(npo_theta, low_eta, donor_alpha=alpha_value)[0]] * len(
            prob_vector)
        y_tf_optimal_dict[alpha_value] = [donor_utility_tf(npo_theta, low_eta, donor_alpha=alpha_value)[1]] * len(
            prob_vector)
        range_dic_l[alpha_value] = [range_dict_l_with_p[k] for k in prob_vector]
        range_dic_h[alpha_value] = [range_dict_h_with_p[k] for k in prob_vector]

    return donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
           npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h


def plots_with_eta(eta_vector, prob_vector=prob, high_eta=eta_h, donor_alpha=alpha,
                   donor_budget=budget,
                   npo_delta=delta, market_beta=beta, npo_theta=theta, low_gamma=gamma_l, num=iter):
    donor_pfr_optimal_dict = {}
    donor_pfr_yt_dict = {}
    u_d_tf_optimal_dict = {}
    y_pfr_optimal_dict = {}
    npo_pfr_optimal_low = {}
    npo_pfr_optimal_high = {}
    y_tf_optimal_dict = {}
    range_dic_l = {}
    range_dic_h = {}
    ratio_h = npo_delta * market_beta / high_eta

    for eta_value in eta_vector:
        ratio_l = npo_delta * market_beta / eta_value
        u_d_pfr_with_p_dict = {}
        yt_with_p_dict = {}
        u_nl_with_p_dict = {}
        u_nh_with_p_dict = {}
        expected_y_with_p_dict = {}
        range_dict_l_with_p = {}
        range_dict_h_with_p = {}
        for p in prob_vector:

            range_dict_l_with_p[p] = ratio_to_region(ratio_l, p, low_gamma)
            range_dict_h_with_p[p] = ratio_to_region(ratio_h, p, low_gamma)
            bar_yt = 2 * max_yt_possible(market_beta, npo_delta, p, low_gamma, high_eta, donor_budget)
            y = np.linspace(0, bar_yt, num)

            u_d_pfr_with_yt_dict = {}
            for j in y:
                u_d_pfr_with_yt_dict[j] = donor_utility_pfr_yt(j, p, low_gamma, eta_value, donor_alpha, npo_theta,
                                                               market_beta, npo_delta,
                                                               donor_budget, high_eta)

            yt_optimal = max(u_d_pfr_with_yt_dict, key=lambda k: u_d_pfr_with_yt_dict[k][0])

            expected_y_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][1]
            u_d_pfr_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][0]  # for a given p, this is optimal u_d_pfr
            u_nl_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][2]
            u_nh_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][3]
            yt_with_p_dict[p] = yt_optimal

        donor_pfr_optimal_dict[eta_value] = [u_d_pfr_with_p_dict[k] for k in prob_vector]
        donor_pfr_yt_dict[eta_value] = [yt_with_p_dict[k] for k in prob_vector]
        y_pfr_optimal_dict[eta_value] = [expected_y_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_low[eta_value] = [u_nl_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_high[eta_value] = [u_nh_with_p_dict[k] for k in prob_vector]
        u_d_tf_optimal_dict[eta_value] = [donor_utility_tf(npo_theta, low_eta=eta_value, donor_alpha=alpha)[0]] * len(
            prob_vector)
        y_tf_optimal_dict[eta_value] = [donor_utility_tf(npo_theta, low_eta=eta_value, donor_alpha=alpha)[1]] * len(
            prob_vector)
        range_dic_l[eta_value] = [range_dict_l_with_p[k] for k in prob_vector]
        range_dic_h[eta_value] = [range_dict_h_with_p[k] for k in prob_vector]

    return donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
           npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h


def plots_with_beta(beta_vector, prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                    donor_alpha=alpha,
                    donor_budget=budget,
                    npo_delta=delta, npo_theta=theta, low_gamma=gamma_l, num=iter):
    donor_pfr_optimal_dict = {}
    donor_pfr_yt_dict = {}
    u_d_tf_optimal_dict = {}
    y_pfr_optimal_dict = {}
    npo_pfr_optimal_low = {}
    npo_pfr_optimal_high = {}
    y_tf_optimal_dict = {}
    range_dic_l = {}
    range_dic_h = {}

    for beta_value in beta_vector:
        ratio_h = npo_delta * beta_value / high_eta
        ratio_l = npo_delta * beta_value / low_eta
        u_d_pfr_with_p_dict = {}
        yt_with_p_dict = {}
        u_nl_with_p_dict = {}
        u_nh_with_p_dict = {}
        expected_y_with_p_dict = {}
        range_dict_l_with_p = {}
        range_dict_h_with_p = {}
        for p in prob_vector:

            range_dict_l_with_p[p] = ratio_to_region(ratio_l, p, low_gamma)
            range_dict_h_with_p[p] = ratio_to_region(ratio_h, p, low_gamma)
            bar_yt = 2 * max_yt_possible(beta_value, npo_delta, p, low_gamma, high_eta, donor_budget)
            y = np.linspace(0, bar_yt, num)

            u_d_pfr_with_yt_dict = {}
            for j in y:
                u_d_pfr_with_yt_dict[j] = donor_utility_pfr_yt(j, p, low_gamma, low_eta, donor_alpha, npo_theta,
                                                               beta_value, npo_delta,
                                                               donor_budget, high_eta)

            yt_optimal = max(u_d_pfr_with_yt_dict, key=lambda k: u_d_pfr_with_yt_dict[k][0])

            expected_y_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][1]
            u_d_pfr_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][0]  # for a given p, this is optimal u_d_pfr
            u_nl_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][2]
            u_nh_with_p_dict[p] = u_d_pfr_with_yt_dict[yt_optimal][3]
            yt_with_p_dict[p] = yt_optimal

        donor_pfr_optimal_dict[beta_value] = [u_d_pfr_with_p_dict[k] for k in prob_vector]
        donor_pfr_yt_dict[beta_value] = [yt_with_p_dict[k] for k in prob_vector]
        y_pfr_optimal_dict[beta_value] = [expected_y_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_low[beta_value] = [u_nl_with_p_dict[k] for k in prob_vector]
        npo_pfr_optimal_high[beta_value] = [u_nh_with_p_dict[k] for k in prob_vector]
        u_d_tf_optimal_dict[beta_value] = [donor_utility_tf(npo_theta, low_eta, donor_alpha)[0]] * len(prob_vector)
        y_tf_optimal_dict[beta_value] = [donor_utility_tf(npo_theta, low_eta, donor_alpha)[1]] * len(prob_vector)
        range_dic_l[beta_value] = [range_dict_l_with_p[k] for k in prob_vector]
        range_dic_h[beta_value] = [range_dict_h_with_p[k] for k in prob_vector]

    return donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
           npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h


def plot_y_lims(donor_pfr_dic, y_pfr_dic, y_tf_dic, npo_pfr_dic, donor_alpha, npo_theta=theta, low_eta=eta_l,
                donor_budget=budget, high_eta=eta_h):
    u_d_max = [donor_utility_tf(npo_theta, high_eta, donor_alpha)[0]]
    u_d_min = [donor_utility_tf(npo_theta, low_eta, donor_alpha)[0]]
    u_n_max = [donor_budget * high_eta]
    u_n_min = -0.25

    y_max = [max(y_tf_dic.values())[0]]
    y_min = [min(y_tf_dic.values())[0]]


    for key in donor_pfr_dic.keys():
        u_d_max.append(max(donor_pfr_dic[key]))
        u_d_min.append(min(donor_pfr_dic[key]))
        u_n_max.append(max(npo_pfr_dic[key]))
        y_max.append(max(y_pfr_dic[key]))
        y_min.append(min(y_pfr_dic[key]))

    u_d_max = max(u_d_max) + 0.25
    u_d_min = min(u_d_min) - 0.25
    y_max = max(y_max) + 0.25
    y_min = min(y_min) - 0.25
    u_n_max = max(u_n_max) + 0.25
    print(y_min)
    return u_d_min, u_d_max, y_min, y_max, u_n_min, u_n_max


def figure_donor_utility(pfr, tf, u_d_min, u_d_max):
    fig, ax = plt.subplots(figsize=[6, 6], facecolor=None, edgecolor='k')
    ax.set_xlabel(r'$ p$', fontsize=30)
    ax.set_ylabel(r'$\bar{u}_d$', fontsize=30)
    ax.scatter(prob, pfr, color='black', marker=".", linewidths=1,
               label=r'$\bar{u}_d^{pfr}$')
    ax.plot(prob, tf, color='black',
            linestyle='solid', linewidth=2, label=r'$\bar{u}_d^{tf}$')
    ax.legend(fontsize=20)
    ax.set_ylim(u_d_min, u_d_max)
    return ax


def figure_y(pfr, tf, y_min, y_max):
    fig, ax = plt.subplots(figsize=[6, 6], facecolor=None, edgecolor='k')
    ax.set_xlabel(r'$ p$', fontsize=30)
    ax.set_ylabel(r'$\bar{y}$', fontsize=30)
    ax.scatter(prob, pfr, color='black', marker=".", linewidths=1,
               label=r'$\bar{y}^{pfr}$')
    ax.plot(prob, tf, color='black', linestyle='solid', linewidth=2,
            label=r'$\bar{y}^{tf}$')
    ax.legend(fontsize=20)
    ax.set_ylim(y_min, y_max)
    return ax


def figure_y_pfr(l, h):
    fig, ax = plt.subplots(figsize=[6, 6], facecolor=None, edgecolor='k')
    ax.set_xlabel(r'$ p$', fontsize=30)
    ax.set_ylabel(r'$\bar{y}$', fontsize=30)
    ax.scatter(prob, l, color='black', marker=".", linewidths=1,
               label=r'$\bar{y}^{l}$')
    ax.plot(prob, h, color='black', linestyle='solid', linewidth=2,
            label=r'$\bar{y}^{h}$')
    ax.legend(fontsize=20)
    return ax


def figure_npo_util(pfr_l, pfr_h, tf_l, tf_h, u_n_min, u_n_max):
    fig, ax = plt.subplots(figsize=[6, 6], facecolor=None, edgecolor='k')
    ax.set_xlabel(r'$ p$', fontsize=30)
    ax.set_ylabel(r'$\bar{u}_n$', fontsize=30)
    #ax.plot(prob, pfr_l, color='black', marker=mark.MarkerStyle(marker='+', fillstyle=None), markevery=5,label=r'$\bar{u}_n{}_{l}^{pfr}$')
    #ax.plot(prob, pfr_h, color='black', marker=mark.MarkerStyle(marker='x', fillstyle=None), markevery=5,label=r'$\bar{u}_n{}_{h}^{pfr}$')
    ax.plot(prob, pfr_l, 'd' , markevery=5,color='black', label=r'$\bar{u}_n{}_{l}^{pfr}$')
    ax.plot(prob, pfr_h, 'x' , markevery=5,color='black', label=r'$\bar{u}_n{}_{h}^{pfr}$')
    ax.plot(prob, tf_l, color='black', linestyle='dashed', linewidth=2,
            label=r'$\bar{u}_n{}_{l}^{tf}$')
    ax.plot(prob, tf_h, color='black', linestyle='dotted', linewidth=2,
            label=r'$\bar{u}_n{}_{h}^{tf}$')
    ax.legend(fontsize=20)
    ax.set_ylim(u_n_min, u_n_max)
    return ax


def plot_graphs_gamma():
    if 2 * eta_bar(theta, eta_l, eta_h) < alpha or (delta * beta < eta_h) \
            or alpha > delta / (math.sqrt(delta * beta / eta_h) - 1):
        print('Error : Parameters out of domain values')
    else:
        donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
        npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h = \
            plots_with_gamma_l(gamma=np.linspace(.25, .75, 3), prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                               donor_budget=budget, npo_delta=delta, market_beta=beta, donor_alpha=alpha,
                               npo_theta=theta)

        npo_tf_optimal_low = [u_n_tf(low_eta=eta_l, donor_alpha=alpha,
                                     npo_theta=theta, high_eta=eta_h, donor_budget=budget)[0]] * len(prob)
        npo_tf_optimal_high = [u_n_tf(low_eta=eta_l, donor_alpha=alpha,
                                      npo_theta=theta, high_eta=eta_h, donor_budget=budget)[1]] * len(prob)

        u_d_min, u_d_max, y_min, y_max, u_n_min, u_n_max = \
            plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high, donor_alpha=alpha,
                        npo_theta=theta, low_eta=eta_l, donor_budget=budget)
        string = ["low", "mid", "high"]

        x = 0
        for key in donor_pfr_optimal_dict.keys():
            print(key, string[x])
            fig = figure_donor_utility(donor_pfr_optimal_dict[key], u_d_tf_optimal_dict[key], u_d_min, u_d_max)

            plt.savefig("u_d_%s_gamma.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y(y_pfr_optimal_dict[key], y_tf_optimal_dict[key], y_min, y_max)
            plt.savefig("y_%s_gamma.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_npo_util(npo_pfr_optimal_low[key], npo_pfr_optimal_high[key], npo_tf_optimal_low,
                                  npo_tf_optimal_high, u_n_min, u_n_max)
            plt.savefig("u_n_%s_gamma.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y_pfr(y_l[key], y_h[key])
            plt.savefig("y_pfr_l_h_%s_gamma.pdf" % string[x], dpi=300, bbox_inches='tight')

            x += 1


def plot_graphs_theta(alpha_local):
    if 2 * eta_bar(theta, eta_l, eta_h) < alpha_local or (delta * beta < eta_h) \
            or alpha_local > delta / (math.sqrt(delta * beta / eta_h) - 1):
        print('Error : Parameters out of domain values')
    else:
        donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
        npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h = \
            plots_with_theta(theta_vector=np.linspace(.25, .75, 3), prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                             donor_budget=budget,
                             npo_delta=delta, market_beta=beta, donor_alpha=alpha_local, low_gamma=gamma_l)

        npo_tf_optimal_low = {}
        npo_tf_optimal_high = {}
        for key in np.linspace(.25, .75, 3):
            npo_tf_optimal_low[key] = [u_n_tf(low_eta=eta_l, donor_alpha=alpha_local,
                                              npo_theta=key, high_eta=eta_h, donor_budget=budget)[0]] * len(prob)
            npo_tf_optimal_high[key] = [u_n_tf(low_eta=eta_l, donor_alpha=alpha_local,
                                               npo_theta=key, high_eta=eta_h, donor_budget=budget)[1]] * len(prob)

        u_d_min, u_d_max, y_min, y_max, u_n_min, u_n_max = \
            plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                    donor_alpha=alpha,
                    npo_theta=theta, low_eta=eta_l, donor_budget=budget)
        if alpha_local < 2 * eta_l:
            alpha_level = "low_alpha"
        else:
            alpha_level = "high_alpha"
        string = ["low", "mid", "high"]

        x = 0
        for key in donor_pfr_optimal_dict.keys():
            fig = figure_donor_utility(donor_pfr_optimal_dict[key], u_d_tf_optimal_dict[key], u_d_min, u_d_max)

            plt.savefig("u_d_%s_theta_%s.pdf" % (string[x], alpha_level), dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y(y_pfr_optimal_dict[key], y_tf_optimal_dict[key], y_min, y_max)
            plt.savefig("y_%s_theta_%s.pdf" % (string[x], alpha_level), dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_npo_util(npo_pfr_optimal_low[key], npo_pfr_optimal_high[key], npo_tf_optimal_low[key],
                                  npo_tf_optimal_high[key], u_n_min, u_n_max)
            plt.savefig("u_n_%s_theta_%s.pdf" % (string[x], alpha_level), dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y_pfr(y_l[key], y_h[key])
            plt.savefig("y_pfr_l_h_%s_theta_%s.pdf" % (string[x], alpha_level), dpi=300, bbox_inches='tight')

            x += 1


def plot_graphs_alpha(alpha_vector):
    if 2 * eta_bar(theta, eta_l, eta_h) < alpha or (delta * beta < eta_h) \
            or alpha > delta / (math.sqrt(delta * beta / eta_h) - 1):
        print('Error : Parameters out of domain values')
    else:
        donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
        npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h, y_l, y_h = \
            plots_with_alpha(alpha_vector, prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                             donor_budget=budget,
                             npo_delta=delta, market_beta=beta, npo_theta=theta, low_gamma=gamma_l)
        print(y_tf_optimal_dict, y_pfr_optimal_dict)
        npo_tf_optimal_low = {}
        npo_tf_optimal_high = {}
        for key in alpha_vector:
            npo_tf_optimal_low[key] = [u_n_tf(low_eta=eta_l, donor_alpha=key,
                                              npo_theta=theta, high_eta=eta_h, donor_budget=budget)[0]] * len(prob)
            npo_tf_optimal_high[key] = [u_n_tf(low_eta=eta_l, donor_alpha=key,
                                               npo_theta=theta, high_eta=eta_h, donor_budget=budget)[1]] * len(prob)

        u_d_min, u_d_max, y_min, y_max, u_n_min, u_n_max = \
            plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                    donor_alpha=alpha,
                    npo_theta=theta, low_eta=eta_l, donor_budget=budget)
        if alpha_vector[0] < 2 * eta_l:
            alpha_level = "low_alpha"
        else:
            alpha_level = "high_alpha"

        string = ["a", "b", "c"]

        x = 0
        for key in donor_pfr_optimal_dict.keys():
            fig = figure_donor_utility(donor_pfr_optimal_dict[key], u_d_tf_optimal_dict[key], u_d_min, u_d_max)

            plt.savefig("figure7%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y(y_pfr_optimal_dict[key], y_tf_optimal_dict[key], y_min, y_max)
            plt.savefig("figure8%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_npo_util(npo_pfr_optimal_low[key], npo_pfr_optimal_high[key], npo_tf_optimal_low[key],
                                  npo_tf_optimal_high[key], u_n_min, u_n_max)
            plt.savefig("figureC10%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            '''
            fig = figure_y_pfr(y_l[key], y_h[key])
            plt.savefig("y_pfr_l_h_%s_alpha_%s.pdf" % (string[x], alpha_level), dpi=300, bbox_inches='tight')
            '''
            x += 1


def plot_graphs_eta(eta_vector=np.linspace(.25, .75, 3)):
    if 2 * eta_bar(theta, eta_vector[0], eta_h) < alpha or (delta * beta < eta_h) \
            or alpha > delta / (math.sqrt(delta * beta / eta_h) - 1):
        print('Error : Parameters out of domain values')
    else:
        donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
        npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h = \
            plots_with_eta(eta_vector, prob_vector=prob, high_eta=eta_h, donor_alpha=alpha,
                           donor_budget=budget,
                           npo_delta=delta, market_beta=beta, npo_theta=theta, low_gamma=gamma_l)

        npo_tf_optimal_low = {}
        npo_tf_optimal_high = {}
        for key in eta_vector:
            npo_tf_optimal_low[key] = [u_n_tf(low_eta=key, donor_alpha=alpha,
                                              npo_theta=theta, high_eta=eta_h, donor_budget=budget)[0]] * len(prob)
            npo_tf_optimal_high[key] = [u_n_tf(low_eta=key, donor_alpha=alpha,
                                               npo_theta=theta, high_eta=eta_h, donor_budget=budget)[1]] * len(prob)

        u_d_min = []
        u_d_max = []
        y_min = []
        y_max = []
        u_n_min = []
        u_n_max = []
        for eta in eta_vector:
            u_d_min.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict,y_tf_optimal_dict, npo_pfr_optimal_high,
                                       donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[0])
            u_d_max.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                                       donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[1])
            y_min.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                                     donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[2])
            y_max.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict,y_tf_optimal_dict, npo_pfr_optimal_high,
                                     donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[3])
            u_n_min.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                                       donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[4])
            u_n_max.append(plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict,y_tf_optimal_dict,  npo_pfr_optimal_high,
                                       donor_alpha=alpha, npo_theta=theta, low_eta=eta, donor_budget=budget)[5])



        u_d_min = min(u_d_min)
        u_d_max = max(u_d_max)
        y_min = min(y_min)
        y_max = max(y_max)
        u_n_min = min(u_n_min)
        u_n_max = max(u_n_max)

        string = ["a", "b", "c"]

        x = 0
        for key in donor_pfr_optimal_dict.keys():
            fig = figure_donor_utility(donor_pfr_optimal_dict[key], u_d_tf_optimal_dict[key], u_d_min, u_d_max)

            plt.savefig("figureB1%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y(y_pfr_optimal_dict[key], y_tf_optimal_dict[key], y_min, y_max)
            plt.savefig("figureB2%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_npo_util(npo_pfr_optimal_low[key], npo_pfr_optimal_high[key], npo_tf_optimal_low[key],
                                  npo_tf_optimal_high[key], u_n_min, u_n_max)
            plt.savefig("figureC11%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            x += 1


def plot_graphs_beta(beta_vector=np.linspace(1.25, 1.75, 3)):
    if 2 * eta_bar(theta, eta_l, eta_h) < alpha or (delta * beta_vector[0] < eta_h) \
            or alpha > delta / (math.sqrt(delta * beta_vector[0] / eta_h) - 1):
        print('Error : Parameters out of domain values')
    else:
        donor_pfr_optimal_dict, u_d_tf_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_low, \
        npo_pfr_optimal_high, donor_pfr_yt_dict, range_dic_l, range_dic_h = \
            plots_with_beta(beta_vector, prob_vector=prob, low_eta=eta_l, high_eta=eta_h,
                            donor_alpha=alpha,
                            donor_budget=budget,
                            npo_delta=delta, npo_theta=theta, low_gamma=gamma_l)

        npo_tf_optimal_low = [u_n_tf(low_eta=eta_l, donor_alpha=alpha,
                                     npo_theta=theta, high_eta=eta_h, donor_budget=budget)[0]] * len(prob)
        npo_tf_optimal_high = [u_n_tf(low_eta=eta_l, donor_alpha=alpha,
                                      npo_theta=theta, high_eta=eta_h, donor_budget=budget)[1]] * len(prob)

        u_d_min, u_d_max, y_min, y_max, u_n_min, u_n_max = \
            plot_y_lims(donor_pfr_optimal_dict, y_pfr_optimal_dict, y_tf_optimal_dict, npo_pfr_optimal_high,
                        donor_alpha=alpha,
                        npo_theta=theta, low_eta=eta_l, donor_budget=budget)
        string = ["a", "b", "c"]

        x = 0
        for key in donor_pfr_optimal_dict.keys():
            print(key, string[x])
            fig = figure_donor_utility(donor_pfr_optimal_dict[key], u_d_tf_optimal_dict[key], u_d_min, u_d_max)

            plt.savefig("figureB7%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_y(y_pfr_optimal_dict[key], y_tf_optimal_dict[key], y_min, y_max)
            plt.savefig("figureB8%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            fig = figure_npo_util(npo_pfr_optimal_low[key], npo_pfr_optimal_high[key], npo_tf_optimal_low,
                                  npo_tf_optimal_high, u_n_min, u_n_max)
            plt.savefig("figureC14%s.pdf" % string[x], dpi=300, bbox_inches='tight')
            plt.show()

            x += 1
