"""Pluggable report composer for GPA and related genomic analysis skills.

Provides a minimal, zero-extra-dependency framework for assembling Markdown
reports from discrete sections. Sections are callables that receive a shared
ReportContext and return either a Markdown string or None (skipped).

Design goals (ponytail full mode):
- One shared context dataclass
- Sections are plain callables / functions
- Composer orders sections, skips empties, joins with newlines
- Built-in registry for common GPA sections
- No template engine beyond Python f-strings
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Protocol, runtime_checkable


@runtime_checkable
class ReportSection(Protocol):
    """A report section callable.

    Implementations receive a ReportContext and return Markdown, or None to be
    omitted from the final report.
    """

    def __call__(self, ctx: "ReportContext") -> Optional[str]: ...


@dataclass
class ReportContext:
    """Shared context passed to every report section.

    Skills may add extra keys via `extras` without changing this dataclass.
    Sections should use getattr for optional fields to stay decoupled.
    """

    variants: List[Any] = field(default_factory=list)
    config: Optional[Any] = None
    tissue_profile: Optional[Dict[str, Any]] = None
    multi_hits: List[Dict[str, Any]] = field(default_factory=list)
    gtex_data: Optional[Dict[str, Dict]] = None
    report_md: Optional[str] = None  # Used by JSON-report sections
    qc_summary: Optional[Dict[str, Any]] = None
    extras: Dict[str, Any] = field(default_factory=dict)

    def get(self, name: str, default: Any = None) -> Any:
        """Fetch an optional attribute or an extra keyed value."""
        if hasattr(self, name):
            return getattr(self, name, default)
        return self.extras.get(name, default)


class ReportComposer:
    """Assemble a Markdown report from ordered ReportSection plugins."""

    def __init__(self, name: str = "report") -> None:
        self.name = name
        self._sections: List[tuple[Optional[float], ReportSection]] = []

    def add_section(
        self,
        section: ReportSection,
        order: Optional[float] = None,
    ) -> "ReportComposer":
        """Register a section. Lower order values appear earlier."""
        self._sections.append((order if order is not None else len(self._sections), section))
        return self

    def add_sections(self, *sections: ReportSection) -> "ReportComposer":
        """Register sections with default ordering."""
        for section in sections:
            self.add_section(section)
        return self

    def compose(self, ctx: ReportContext, joiner: str = "\n") -> str:
        """Run all sections, skip empties, and join the Markdown output."""
        ordered = sorted(self._sections, key=lambda x: x[0])
        parts: List[str] = []
        for _, section in ordered:
            try:
                text = section(ctx)
            except Exception as exc:
                # A section failure must not kill the whole report.
                text = f"<!-- section error: {type(exc).__name__}: {exc} -->"
            if text:
                parts.append(text.rstrip())
        return joiner.join(parts)

    @property
    def section_count(self) -> int:
        return len(self._sections)


# ---------------------------------------------------------------------------
# Built-in generic sections that any skill can reuse
# ---------------------------------------------------------------------------


def header_section(title: str, subtitle: Optional[str] = None) -> ReportSection:
    """Return a section that emits a simple markdown header block."""

    def _section(ctx: ReportContext) -> str:
        lines = [f"# {title}"]
        if subtitle:
            lines.append(f"\n{subtitle}")
        return "\n".join(lines)

    return _section


def static_section(text: str) -> ReportSection:
    """Return a section that always emits the given static markdown."""

    def _section(ctx: ReportContext) -> str:
        return text

    return _section
