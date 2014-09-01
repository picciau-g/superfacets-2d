/****************************************************************************
** Meta object code from reading C++ file 'dialogs.h'
**
** Created: Wed Apr 16 00:18:51 2014
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "dialogs.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'dialogs.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_dialogS[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
       9,    8,    8,    8, 0x0a,
      21,    8,    8,    8, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_dialogS[] = {
    "dialogS\0\0nameGiven()\0buttonSaveClicked()\0"
};

void dialogS::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        dialogS *_t = static_cast<dialogS *>(_o);
        switch (_id) {
        case 0: _t->nameGiven(); break;
        case 1: _t->buttonSaveClicked(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData dialogS::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject dialogS::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_dialogS,
      qt_meta_data_dialogS, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &dialogS::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *dialogS::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *dialogS::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_dialogS))
        return static_cast<void*>(const_cast< dialogS*>(this));
    return QWidget::qt_metacast(_clname);
}

int dialogS::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
