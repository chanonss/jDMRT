'use strict';

// Same model-terms / post hoc behaviour as jamovi's own ANOVA analysis.
module.exports = {

    view_updated: function(ui) {
        this.calcModelTerms(ui);
        this.filterModelTerms(ui);
    },

    factors_changed: function(ui) {
        this.calcModelTerms(ui);
    },

    blk_changed: function(ui) {
        this.calcModelTerms(ui);
    },

    // the block variable (if any) is offered as a model component too, but is never crossed
    // with the fixed factors automatically
    blockName: function(ui) {
        let b = ui.blk ? ui.blk.value() : null;
        if (Array.isArray(b))
            b = b.length > 0 ? b[0] : null;
        return (typeof b === 'string' && b.length > 0) ? b : null;
    },

    modelTerms_changed: function(ui) {
        this.filterModelTerms(ui);
    },

    modelSupplier_updated: function(ui) {
        let variableList = utils.clone(ui.factors.value(), []);
        let blk = this.blockName(ui);
        if (blk !== null)
            variableList.push(blk);
        ui.modelSupplier.setValue(utils.valuesToItems(variableList, FormatDef.variable));
    },

    postHocSupplier_updated: function(ui) {
        let termsList = utils.clone(ui.modelTerms.value(), []);
        ui.postHocSupplier.setValue(utils.valuesToItems(termsList, FormatDef.term));
    },

    postHocSupplier_changed: function(ui) {
        let values = utils.itemsToValues(ui.postHocSupplier.value());
        utils.checkValue(ui.postHoc, true, values, FormatDef.term);
    },

    filterModelTerms: function(ui) {
        var termsList = utils.clone(ui.modelTerms.value(), []);

        // remove terms that contain a term which was just removed
        var termsDiff = this.findChanges("currentList", termsList, true, FormatDef.term);
        var changed = false;
        if (termsDiff.removed.length > 0 && termsList !== null) {
            var itemsRemoved = false;
            for (var i = 0; i < termsDiff.removed.length; i++) {
                var item = termsDiff.removed[i];
                for (var j = 0; j < termsList.length; j++) {
                    if (FormatDef.term.contains(termsList[j], item)) {
                        termsList.splice(j, 1);
                        j -= 1;
                        itemsRemoved = true;
                    }
                }
            }
            if (itemsRemoved)
                changed = true;
        }

        if (utils.sortArraysByLength(termsList))
            changed = true;

        if (changed)
            ui.modelTerms.setValue(termsList);

        ui.postHocSupplier.setValue(utils.valuesToItems(termsList, FormatDef.term));
    },

    calcModelTerms: function(ui) {
        var variableList = utils.clone(ui.factors.value(), []);
        var blk = this.blockName(ui);
        var supply = variableList.slice();
        if (blk !== null)
            supply.push(blk);

        ui.modelSupplier.setValue(utils.valuesToItems(supply, FormatDef.variable));

        var varsDiff = this.findChanges("variableList", variableList, true, FormatDef.variable);
        var termsList = utils.clone(ui.modelTerms.value(), []);

        var termsChanged = false;
        for (var i = 0; i < varsDiff.removed.length; i++) {
            for (var j = 0; j < termsList.length; j++) {
                if (FormatDef.term.contains(termsList[j], varsDiff.removed[i])) {
                    termsList.splice(j, 1);
                    termsChanged = true;
                    j -= 1;
                }
            }
        }

        termsList = utils.getCombinations(varsDiff.added, termsList);
        termsChanged = termsChanged || varsDiff.added.length > 0;

        // drop terms that use a variable which is neither a factor nor the current block
        for (var j = 0; j < termsList.length; j++) {
            var stale = false;
            for (var k = 0; k < termsList[j].length; k++)
                if (supply.indexOf(termsList[j][k]) === -1)
                    stale = true;
            if (stale) {
                termsList.splice(j, 1);
                j -= 1;
                termsChanged = true;
            }
        }

        // the block enters the model as one additive term
        if (blk !== null) {
            var has = false;
            for (var m = 0; m < termsList.length; m++)
                if (termsList[m].length === 1 && termsList[m][0] === blk)
                    has = true;
            if (!has) {
                termsList.push([blk]);
                termsChanged = true;
            }
        }

        if (termsChanged)
            ui.modelTerms.setValue(termsList);
    }

};
