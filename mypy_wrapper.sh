#!/usr/bin/env bash
/opt/local/Library/Frameworks/Python.framework/Versions/3.6/bin/mypy "$@" |  awk -v prefix="/" '{print prefix $0}'

