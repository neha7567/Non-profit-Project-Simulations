import math
import numpy as np

budget = 1
eta_h = 1
eta_l = 0.5
theta = 0.50
delta = 2
beta = 1.2
alpha = .75
gamma_l = 0.50
prob = np.linspace(.01, .99, 100)


def higher_shock(p, low_gamma):
    try:
        gamma_h = (1 - p * low_gamma) / (1 - p)
    except ZeroDivisionError:
        gamma_h = 1
    return gamma_h


# Given a target y_T, the maximum npo's utility and corresponding y delivered.
def find_expected_y_l(yt, p, low_gamma, market_beta, low_eta, npo_delta=delta, donor_budget=budget):
    gamma_h = higher_shock(p, low_gamma)
    region_utility_dict = {}
    region_y_dict = {}
    ratio = npo_delta * market_beta / low_eta