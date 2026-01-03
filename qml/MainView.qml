import QtQuick 6.4
import QtQuick.Controls 6.4

ApplicationWindow {
    width: 640
    height: 480
    visible: true
    color: "black"

    property int rectCount: 24

    Row {
        id: row
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 16
        height: 80
        spacing: 8

        Repeater {
            model: rectCount
            delegate: Rectangle {
                width: Math.max(20, (row.width - (rectCount - 1) * row.spacing) / rectCount)
                height: row.height
                color: "#e74c3c"
                radius: 6

                MouseArea {
                    anchors.fill: parent
                    onClicked: ui.rectangleClicked(index)
                }
            }
        }
    }
}
