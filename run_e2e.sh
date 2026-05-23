#!/bin/bash
cd /root/.openclaw/skills/dgra-genomic-risk
export PYTHONPATH=scripts
python3 tests/test_v0_9_3_e2e.py
