install-me:
	swipl -f none -g "pack_install(.,[upgrade(true)]), halt"

publish:
	swipl -f none -g "pack_property(prism,download(D)), pack_install(D), halt"

all:
check:
install:
