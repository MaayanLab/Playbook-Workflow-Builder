## Playbook Partnership Installation Guide
When starting from scratch (without nodejs, python, git, or vscode installed), this document may help you get started with the prerequisites to the [contributions guide](./contributions.md) depending on your operating system. It is not the only way these things can be installed but it's a way which has been tested to work.

### Windows

#### With Scoop Installed

See <https://scoop.sh/> for install instructions.

From a powershell prompt run:
```powershell
scoop bucket add main
scoop install nodejs git python
scoop bucket add extras
scoop install vscode
```

#### With Chocolatey Insalled

See <https://chocolatey.org/install> for install instructions.

From an admin powershell prompt run:
```powershell
choco install vscode nodejs git python==3.8.5 -y
```

### Mac OS X
1. Install Homebrew: https://brew.sh/
2. From console run:
```
brew install git node python@3.8
brew install --cask visual-studio-code
```

### Ubuntu
1. From console run:
```bash
sudo apt update
sudo apt install git nodejs
sudo snap install --classic code
```
2. Setup local npm install, add to PATH. This is an opinionated setup which stores all npm globally installed packages (`npm i -g`) in your local profile (avoiding permission issues). The `~/.profile` is assumed to be sourced at startup by your shell, depending on your system it may need to be renamed to something else like `~/.bashrc` or `~/.zshrc`.
```bash
mkdir -p $HOME/.npm-global
npm config set prefix=$HOME/.npm-global
npm i -g npm
echo "PATH=$HOME/.npm-global/bin:$PATH" >> ~/.profile
source ~/.profile
npm i -g node
```
