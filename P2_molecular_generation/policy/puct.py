"""
PUCT (Predictor + Upper Confidence bound for Trees) policy.

AlphaGo-style selection that incorporates RNN prior probability
into the UCB formula, guiding early exploration toward tokens
the RNN considers promising.

PUCT = Q(s,a)/N(s,a) + c_puct * P(s,a) * sqrt(N(s)) / (1 + N(s,a))

where P(s,a) is the RNN's predicted probability for token a given state s.
"""

from math import sqrt

from chemtsv2.abc import Policy


class Puct(Policy):
    def evaluate(child_state, conf):
        c_puct = conf.get("c_puct", 1.0)

        # Get prior probability stored on the child node
        prior = getattr(child_state, "prior", 0.0)

        exploitation = child_state.total_reward / child_state.visits
        exploration = c_puct * prior * sqrt(child_state.parent_node.state.visits) / (1 + child_state.visits)

        return exploitation + exploration
