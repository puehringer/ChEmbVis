FROM node:12-alpine
RUN mkdir /frontend
WORKDIR /frontend
# Install app dependencies
COPY package.json ./
COPY package-lock.json ./
RUN npm install
# Add app
COPY . ./
# Run the dev server
CMD ["npm", "start"]
