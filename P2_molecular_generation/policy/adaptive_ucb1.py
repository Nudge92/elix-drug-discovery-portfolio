"""
Adaptive UCB1 policy: c_val decays from c_start to c_end over the course of generation.

Early phase: high c_val (exploration) → discover diverse chemical space
Late phase:  low c_val (exploitation) → refine best regions found

Config parameters:
  c_start: initial c_val (default: 1.0)
  c_end:   final c_val (default: 0.1)
  adaptive_schedule: "linear" or "step" (default: "linear")
    - linear: c = c_start + (c_end - c_start) * progress
    - step:   c = c_start if progress < 0.5 else c_end
"""

from math import log, sqrt

from chemtsv2.abc import Policy


class AdaptiveUcb1(Policy):
    def evaluate(child_state, conf):
        c_start = conf.get("c_start", 1.0)
        c_end = conf.get("c_end", 0.1)
        schedule = conf.get("adaptive_schedule", "linear")

        # Progress: fraction of generation completed (0.0 → 1.0)
        progress = conf.get("_progress", 0.0)

        if schedule == "linear":
            c_val = c_start + (c_end - c_start) * progress
        elif schedule == "step":
            c_val = c_start if progress < 0.5 else c_end
        else:
            c_val = c_start

        ucb1 = (child_state.total_reward / child_state.visits) + c_val * sqrt(
            2 * log(child_state.parent_node.state.visits) / child_state.visits
        )
        return ucb1
