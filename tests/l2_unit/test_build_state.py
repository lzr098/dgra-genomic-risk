#!/usr/bin/env python3
"""L2 unit tests for dgra_build_state.py"""
import pytest
from pathlib import Path
import json

import dgra_build_state as bs


class TestLoadSaveState:
    def test_load_missing_returns_empty(self, tmp_path):
        bs.set_state_file(tmp_path / "nonexistent.json")
        assert bs.load_state() == {}

    def test_save_and_load(self, tmp_path):
        state_file = tmp_path / "state.json"
        bs.set_state_file(state_file)
        bs.save_state("test_step", {"status": "complete", "foo": "bar"})
        state = bs.load_state()
        assert "test_step" in state
        assert state["test_step"]["status"] == "complete"
        assert state["test_step"]["data"]["foo"] == "bar"

    def test_load_corrupted_json_returns_empty(self, tmp_path):
        state_file = tmp_path / "state.json"
        state_file.write_text("not json")
        bs.set_state_file(state_file)
        assert bs.load_state() == {}

    def test_multiple_steps(self, tmp_path):
        state_file = tmp_path / "state.json"
        bs.set_state_file(state_file)
        bs.save_state("step1", {"status": "complete"})
        bs.save_state("step2", {"status": "in_progress"})
        state = bs.load_state()
        assert state["step1"]["status"] == "complete"
        assert state["step2"]["status"] == "in_progress"


class TestGetStepStatus:
    def test_pending_when_missing(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        assert bs.get_step_status("missing") == "pending"

    def test_returns_status(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step", {"status": "complete"})
        assert bs.get_step_status("step") == "complete"


class TestGetStepData:
    def test_none_when_missing(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        assert bs.get_step_data("missing") is None

    def test_returns_data(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step", {"status": "complete", "extra": 123})
        assert bs.get_step_data("step")["extra"] == 123


class TestIsStepComplete:
    def test_true_when_complete(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step", {"status": "complete"})
        assert bs.is_step_complete("step") is True

    def test_false_when_in_progress(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step", {"status": "in_progress"})
        assert bs.is_step_complete("step") is False


class TestResetState:
    def test_reset_single_step(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step1", {"status": "complete"})
        bs.save_state("step2", {"status": "complete"})
        bs.reset_state("step1")
        state = bs.load_state()
        assert "step1" not in state
        assert "step2" in state

    def test_reset_all_deletes_file(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step", {"status": "complete"})
        bs.reset_state()
        assert not (tmp_path / "state.json").exists()


class TestListCompletedSteps:
    def test_lists_only_complete(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        bs.save_state("step1", {"status": "complete"})
        bs.save_state("step2", {"status": "failed"})
        bs.save_state("step3", {"status": "complete"})
        completed = bs.list_completed_steps()
        assert sorted(completed) == ["step1", "step3"]


class TestBuildStepContextManager:
    def test_successful_completion(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        with bs.BuildStep("my_step") as step:
            step.complete(result="ok")
        assert bs.is_step_complete("my_step")
        assert bs.get_step_data("my_step")["result"] == "ok"

    def test_exception_records_failure(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        with pytest.raises(ValueError):
            with bs.BuildStep("fail_step") as step:
                raise ValueError("boom")
        assert bs.get_step_status("fail_step") == "failed"

    def test_in_progress_on_enter(self, tmp_path):
        bs.set_state_file(tmp_path / "state.json")
        with bs.BuildStep("running_step") as step:
            assert bs.get_step_status("running_step") == "in_progress"
            step.complete()
