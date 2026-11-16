#!/bin/bash
cd /home/ubuntu/app
exec npx tsx cli/cavatica-worker.ts "$@"
