import gntp.notifier

image = open('img.png').read()
growl = gntp.notifier.GrowlNotifier(
    applicationName="Image downloader",
    notifications=["New Messages"],
    defaultNotifications=["New Messages"],
)
growl.register()


def notify(title, description):
    growl.notify(
        noteType="New Messages",
        title=title,
        description=description,
        icon=image,
        sticky=False,
        priority=1,
    )
