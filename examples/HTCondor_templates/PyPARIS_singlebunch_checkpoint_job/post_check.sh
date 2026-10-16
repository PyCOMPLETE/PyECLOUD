#!/bin/bash
ec="$1"

case "$1" in
  0)   exit 0 ;;    # finished
  177) exit 177 ;;  # retry
  *)   exit 11 ;;   # fatal error
esac