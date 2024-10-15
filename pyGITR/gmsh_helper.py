# -*- coding: utf-8 -*-


def SetGroups(Model, Tags, Name, Color=None):
    def CleanEntities(Model, Tag, GenerateListEntities = True):
        LTag = []
        for t in Tag:
            if type(t) == list:
                LTag.extend(t)
            else:
                LTag.append(t)
        Tag = LTag
        ListEntities = Model.getEntities(2)
        Tag = [(2,t) if type(t) == int else t for t in Tag]
        if GenerateListEntities:
            return [t[1] for t in Tag if t in ListEntities]
        else:
            return [t for t in Tag if t in ListEntities]

    Group = Model.addPhysicalGroup(dim=2, tags=CleanEntities(Model, Tags))
    Model.setPhysicalName(dim=2, tag=Group, name=Name)
    Tag = [(2,g) for g in Model.getEntitiesForPhysicalGroup(2,Group)]
    if Color is not None:
        Model.setColor(Tag, Color[0], Color[1], Color[2])