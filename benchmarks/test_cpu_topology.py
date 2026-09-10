from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from benchmarks import cpu_topology


class CpuTopologyTests(unittest.TestCase):
    def create_cpu(
        self,
        sysfs: Path,
        cpu: int,
        core: int,
        *,
        package: int = 0,
        die: int | None = None,
    ) -> None:
        topology = sysfs / f"cpu{cpu}" / "topology"
        topology.mkdir(parents=True)
        (topology / "physical_package_id").write_text(str(package), encoding="utf-8")
        (topology / "core_id").write_text(str(core), encoding="utf-8")
        if die is not None:
            (topology / "die_id").write_text(str(die), encoding="utf-8")

    def detect(
        self,
        sysfs: Path,
        allowed: set[int],
        *,
        slurm_limit: str | None = None,
    ) -> cpu_topology.CpuTopology:
        environment = (
            {} if slurm_limit is None else {"SLURM_CPUS_PER_TASK": slurm_limit}
        )
        with (
            mock.patch.dict(os.environ, environment, clear=True),
            mock.patch.object(
                os, "sched_getaffinity", return_value=allowed, create=True
            ),
            mock.patch.object(cpu_topology, "SYSFS_CPU_ROOT", sysfs),
            mock.patch.object(cpu_topology, "cpu_model", return_value="Test processor"),
        ):
            return cpu_topology.detect_cpu_topology()

    def test_hybrid_20_core_28_cpu_order_puts_physical_cores_before_smt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            sysfs = Path(directory)
            for cpu in range(28):
                self.create_cpu(sysfs, cpu, cpu // 2 if cpu < 16 else cpu - 8)
            result = self.detect(sysfs, set(range(28)))

        self.assertEqual(result.logical_cpus, 28)
        self.assertEqual(result.physical_cores, 20)
        self.assertEqual(
            result.cpu_order,
            (*range(0, 16, 2), *range(16, 28), *range(1, 16, 2)),
        )
        self.assertEqual(result.affinity_cpus, tuple(range(28)))
        self.assertEqual(result.model, "Test processor")
        self.assertIsNone(result.slurm_cpus_per_task)
        self.assertIn("sysfs", result.source)

    def test_sparse_affinity_counts_only_permitted_cores_and_siblings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            sysfs = Path(directory)
            for cpu in range(8):
                self.create_cpu(sysfs, cpu, cpu % 4)
            result = self.detect(sysfs, {1, 3, 4, 5, 7})

        self.assertEqual(result.logical_cpus, 5)
        self.assertEqual(result.physical_cores, 3)
        self.assertEqual(result.affinity_cpus, (1, 3, 4, 5, 7))
        self.assertEqual(result.cpu_order, (1, 3, 4, 5, 7))

    def test_identical_core_ids_on_distinct_packages_and_dies_stay_distinct(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            sysfs = Path(directory)
            for cpu, package, die in ((0, 0, 0), (1, 0, 0), (2, 1, 0), (3, 0, 1)):
                self.create_cpu(sysfs, cpu, 0, package=package, die=die)
            result = self.detect(sysfs, {0, 1, 2, 3})

        self.assertEqual(result.physical_cores, 3)
        self.assertEqual(result.cpu_order, (0, 2, 3, 1))

    def test_incomplete_topology_preserves_sorted_affinity_without_core_estimate(
        self,
    ) -> None:
        for invalid in (None, "not an integer", "-1"):
            with (
                self.subTest(core_id=invalid),
                tempfile.TemporaryDirectory() as directory,
            ):
                sysfs = Path(directory)
                self.create_cpu(sysfs, 1, 0)
                if invalid is not None:
                    self.create_cpu(sysfs, 5, 1)
                    (sysfs / "cpu5/topology/core_id").write_text(
                        invalid, encoding="utf-8"
                    )
                result = self.detect(sysfs, {5, 1})
            self.assertEqual(result.logical_cpus, 2)
            self.assertIsNone(result.physical_cores)
            self.assertEqual(result.cpu_order, (1, 5))
            self.assertIn("unavailable", result.source)

    def test_slurm_cap_selects_physical_cores_before_siblings(self) -> None:
        for limit, expected_order, physical in (
            ("2", (0, 2), 2),
            ("3", (0, 2, 1), 2),
            ("8", (0, 2, 1, 3), 2),
        ):
            with (
                self.subTest(limit=limit),
                tempfile.TemporaryDirectory() as directory,
            ):
                sysfs = Path(directory)
                for cpu in range(4):
                    self.create_cpu(sysfs, cpu, cpu // 2)
                result = self.detect(sysfs, {0, 1, 2, 3}, slurm_limit=limit)
            self.assertEqual(result.logical_cpus, len(expected_order))
            self.assertEqual(result.cpu_order, expected_order)
            self.assertEqual(result.physical_cores, physical)
            self.assertEqual(result.affinity_cpus, (0, 1, 2, 3))
            self.assertEqual(result.slurm_cpus_per_task, int(limit))
            self.assertIn("SLURM_CPUS_PER_TASK", result.source)

    def test_partially_available_die_ids_do_not_overcount_physical_cores(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            sysfs = Path(directory)
            self.create_cpu(sysfs, 0, 0, die=0)
            self.create_cpu(sysfs, 1, 0)
            result = self.detect(sysfs, {0, 1})

        self.assertIsNone(result.physical_cores)
        self.assertEqual(result.cpu_order, (0, 1))

    def test_malformed_slurm_limits_are_actionable_errors(self) -> None:
        for limit in ("", "0", "-1", "1.5", "4(x2)", "many", " 2", "²"):
            with (
                self.subTest(limit=limit),
                mock.patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": limit}, clear=True),
                self.assertRaisesRegex(
                    ValueError, "SLURM_CPUS_PER_TASK.*positive integer"
                ),
            ):
                cpu_topology.detect_cpu_topology()

    def test_empty_affinity_is_an_error(self) -> None:
        with (
            tempfile.TemporaryDirectory() as directory,
            self.assertRaisesRegex(ValueError, "affinity.*no available CPUs"),
        ):
            self.detect(Path(directory), set())

    def test_fallback_prefers_process_count_and_applies_scheduler_limit(self) -> None:
        with (
            mock.patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "3"}, clear=True),
            mock.patch.object(
                os, "sched_getaffinity", side_effect=OSError, create=True
            ),
            mock.patch.object(os, "process_cpu_count", return_value=6, create=True),
            mock.patch.object(os, "cpu_count", return_value=12) as host_count,
            mock.patch.object(cpu_topology, "cpu_model", return_value="Portable CPU"),
        ):
            result = cpu_topology.detect_cpu_topology()

        self.assertEqual(result.logical_cpus, 3)
        self.assertIsNone(result.physical_cores)
        self.assertIsNone(result.cpu_order)
        self.assertIsNone(result.affinity_cpus)
        self.assertIn("os.process_cpu_count", result.source)
        host_count.assert_not_called()

    def test_fallback_uses_host_count_when_process_count_is_unavailable(self) -> None:
        for process_count in (None, mock.Mock(return_value=None)):
            with (
                self.subTest(process_count=process_count),
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(os, "sched_getaffinity", None, create=True),
                mock.patch.object(os, "process_cpu_count", process_count, create=True),
                mock.patch.object(os, "cpu_count", return_value=8),
                mock.patch.object(
                    cpu_topology, "cpu_model", return_value="Portable CPU"
                ),
            ):
                result = cpu_topology.detect_cpu_topology()
            self.assertEqual(result.logical_cpus, 8)
            self.assertIsNone(result.cpu_order)
            self.assertIsNone(result.physical_cores)
            self.assertEqual(result.source, "os.cpu_count")

    def test_unavailable_or_empty_fallback_counts_are_errors(self) -> None:
        for count in (None, 0, -1):
            with (
                self.subTest(count=count),
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(os, "sched_getaffinity", None, create=True),
                mock.patch.object(os, "process_cpu_count", None, create=True),
                mock.patch.object(os, "cpu_count", return_value=count),
                self.assertRaisesRegex(ValueError, "available CPU count"),
            ):
                cpu_topology.detect_cpu_topology()


if __name__ == "__main__":
    unittest.main()
