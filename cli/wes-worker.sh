#!/bin/bash
cd /home/ubuntu/app
npm start -- "{\"port\":3001,\"plugins\":[\"next\",\"ws\",\"cavatica-proxy\"],\"proxy\":$@}"
