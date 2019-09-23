pyinstaller --noconfirm --log-level=WARN \
    --onefile --nowindow \
    --hidden-import=traits \
    --hidden-import=traitsui \
    --upx-dir=/usr/local/share/ \
    main.spec
