"""K1 回归测试：三个预测/迁移基线的 immigrant_ratio 必须真实生效。

修复前，AR1 / GMM / PPS 按「平移/采样恰好 N 行 + 末尾 append 移民 +
_resize 截断前 N 行」的顺序运行，追加的随机移民被全部截断丢弃
（实际行为为纯预测迁移、无移民）。本组测试锁定修复后的行为：
新种群 = 前 N-n_imm 行平移/采样结果 + n_imm 行均匀随机移民。

测试原理：用同一随机种子分别构造 ratio=0（纯预测）与默认 ratio
两个实例；修复正确时两者前 N-n_imm 行完全一致、移民区必然不同。
若回退到旧实现（移民全丢），则两个输出完全相同，断言失败。
"""
import numpy as np
import pytest

from algorithms import ARNSGA2, GMMTransfer, PPSNSGA2
from problems import FDA1


def _make_case(n=30):
    p = FDA1(tau_t=5, n_t=10)
    p.set_generation(5)
    rng = np.random.default_rng(123)
    X = rng.uniform(p.xl, p.xu, size=(n, p.n_var))
    F = p._calc_objectives(X)
    return p, X, F


@pytest.mark.parametrize("cls", [ARNSGA2, PPSNSGA2],
                         ids=["ar1_predictive", "pps"])
def test_prediction_baseline_keeps_immigrants(cls):
    p, X, F = _make_case()
    n = X.shape[0]
    n_imm = int(round(0.2 * n))
    X0, _ = cls(immigrant_ratio=0.0).respond(
        X, F, p, 1, np.random.default_rng(0))
    X2, _ = cls(immigrant_ratio=0.2).respond(
        X, F, p, 1, np.random.default_rng(0))
    assert X0.shape == (n, p.n_var) and X2.shape == (n, p.n_var)
    # 平移部分保留在原位（与纯预测结果的前 N-n_imm 行一致）
    assert np.array_equal(X2[:n - n_imm], X0[:n - n_imm])
    # 移民区为真实随机移民（旧 bug 下 X2==X0，此断言必失败）
    assert not np.array_equal(X2[n - n_imm:], X0[n - n_imm:])


def test_gmm_baseline_keeps_immigrants():
    p, X, F = _make_case()
    n = X.shape[0]
    n_imm = int(round(0.1 * n))
    X0, _ = GMMTransfer(immigrant_ratio=0.0).respond(
        X, F, p, 1, np.random.default_rng(0))
    X2, _ = GMMTransfer(immigrant_ratio=0.1).respond(
        X, F, p, 1, np.random.default_rng(0))
    assert X0.shape == (n, p.n_var) and X2.shape == (n, p.n_var)
    # GMM 采样部分保留（gmm_seed 消耗 rng 的方式两实例一致）
    assert np.array_equal(X2[:n - n_imm], X0[:n - n_imm])
    # 移民区为真实随机移民（旧 bug 下 X2==X0，此断言必失败）
    assert not np.array_equal(X2[n - n_imm:], X0[n - n_imm:])
