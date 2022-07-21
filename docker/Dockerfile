FROM debian:stable-slim

RUN apt-get update -y && apt-get -y install sudo git make gcc autoconf pkg-config xorg-dev x11-apps exuberant-ctags imagemagick python3 python3-matplotlib python3-pandas r-base gnuplot tzdata python-is-python3 zsh texlive-latex-base dvipng texlive-latex-extra texlive-fonts-recommended cm-super

RUN chsh -s $(which zsh)
RUN echo "alias ls='ls -FC --color=auto'" > $HOME/.zshrc \
      && echo "alias ll='ls -lg'" >> $HOME/.zshrc \
      && echo "alias la='ls -A'"  >> $HOME/.zshrc \
      && echo "alias -g G='| grep --color=auto'" >> $HOME/.zshrc \
      && echo "autoload history-search-end" >> $HOME/.zshrc \
      && echo "zle -N history-beginning-search-backward-end history-search-end" >> $HOME/.zshrc \
      && echo "zle -N history-beginning-search-forward-end history-search-end" >> $HOME/.zshrc \
      && echo 'bindkey "^P" history-beginning-search-backward-end' >> $HOME/.zshrc \
      && echo 'bindkey "^N" history-beginning-search-forward-end' >> $HOME/.zshrc

WORKDIR $HOME
CMD ["zsh"]
