import QtQuick 6.4
import QtQuick.Controls 6.4

ApplicationWindow {
    width: 640
    height: 480
    visible: true
    color: "black"

    property int rectCount: Math.max(1, ui.binCount)
    property int minRectWidth: 8
    property int activeIndex: -1

    Item {
        id: binStrip
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 16
        height: 80

        Row {
            id: row
            anchors.fill: parent
            spacing: 4

            property real cellWidth: Math.max(minRectWidth, (row.width - (rectCount - 1) * row.spacing) / rectCount)

            Repeater {
                model: rectCount
                delegate: Rectangle {
                    width: row.cellWidth
                    height: row.height
                    color: index === activeIndex ? "#f39c12" : "#e74c3c"
                    radius: 6
                }
            }
        }

        MouseArea {
            anchors.fill: parent
            hoverEnabled: false

            function updateIndex(xPos) {
                var cell = row.cellWidth + row.spacing
                var idx = Math.floor(xPos / cell)
                if (idx < 0 || idx >= rectCount) {
                    return
                }
                if (activeIndex !== idx) {
                    activeIndex = idx
                    ui.rectangleClicked(idx)
                }
            }

            onPressed: updateIndex(mouse.x)
            onPositionChanged: if (pressed) updateIndex(mouse.x)
        }
    }
}
