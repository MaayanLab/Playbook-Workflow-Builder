## PREPARE: This step gathers only dependency specification files for
#           ideal caching properties (only reinstall if deps change)
FROM node:latest as prepare

WORKDIR /app

COPY . .
RUN set -x \
  find /app \
    \! \( -name "package.json" -o -name "requirements.txt" \) \
  | xargs rm -rf

FROM node:latest

WORKDIR /app
EXPOSE 3000

# we install necessary system dependencies
RUN set -x \
  && apt-get -y update \
  && apt-get -y install \
    python3-dev \
    python3-pip

# we bring in prepared depedency files
COPY --from=prepare /app .
RUN set -x \
  && npm i \
  && pip install -r requirements.txt

# we copy the rest over to build
ADD . /app
RUN npm run build

CMD ["npm", "run", "start"]