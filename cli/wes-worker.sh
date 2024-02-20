#!/bin/sh
cd /app
npm start -- "{\"port\":3001,\"plugins\":[\"next\",\"ws\",\"cavatica-proxy\"],\"proxy\":$@}"
