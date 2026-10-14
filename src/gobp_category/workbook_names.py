"""Excel-safe names for Group and comparison outputs."""

import re


def unique_sheet_title(value: str, used: set) -> str:
    """Return a distinct Excel worksheet title (31 characters maximum)."""
    base = re.sub(r"[\\/*?:\[\]]", "_", str(value)).strip(" '") or "Group"
    base = base[:31]
    title = base
    index = 2
    while title.casefold() in used:
        suffix = f"_{index}"
        title = f"{base[:31 - len(suffix)]}{suffix}"
        index += 1
    used.add(title.casefold())
    return title


def unique_group_filename(value: str, used: set) -> str:
    """Return a portable filename component while retaining the Group label."""
    base = re.sub(r'[\\/:*?"<>|\x00-\x1f]', "_", str(value)).strip(" .") or "Group"
    base = base[:80]
    name = base
    index = 2
    while name.casefold() in used:
        suffix = f"_{index}"
        name = f"{base[:80 - len(suffix)]}{suffix}"
        index += 1
    used.add(name.casefold())
    return name
