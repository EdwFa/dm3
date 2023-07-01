import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';

import { variables } from './Variables.js';


export class Main extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.state = {
            updateOr: false,
            loading: false,
            useAll: true,

            articles: [],
            token: variables.token,
            DetailArticle: null,
            articlesInfo: [],
            translation_stack: null,
            full_query: null,
            short_query: null,
            task: null,
            message: null,

            // Filters
            queryText: null,
            queryStartDate: null,
            queryEndDate: null,
            queryTypes: new Set(),
            queryOlds: new Set(),
            queryGenders: new Set(),

            // Graph
            graph: null,
        }
    }

    createQuery() {
        let query = "/"
        if (this.state.queryText) {
            query = `${query}?search_field=${this.state.queryText}`
        }
        if (this.state.queryStartDate) {
            query = `${query}&dateStart=${this.state.queryStartDate}`
        }
        if (this.state.queryEndDate) {
            query = `${query}&dateStop=${this.state.queryEndDate}`
        }
        if (this.state.queryGenders.size != 0) {
            for (let el of [...this.state.queryGenders]) {
                query = `${query}&Gender=${el}`
            }
        }
        if (this.state.queryTypes.size != 0) {
            for (let el of [...this.state.queryTypes]) {
                query = `${query}&Type=${el}`
            }
        }
        if (this.state.queryOlds.size != 0) {
            for (let el of [...this.state.queryOlds]) {
                query = `${query}&Old=${el}`
            }
        }

        if (query !== "/" & !this.state.queryText) {
            throw "Not query text"
        }

        return query
    }

    getArticles = (url, interval = 10000) => {
      fetch(variables.API_URL + '/articles',
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${variables.token}`,
                },
            }
          )
          .then(response => {
                console.log(response.status);
                if (response.status == 403) {
                    setTimeout(() => {
                      return this.getArticles(url, interval)
                    }, interval);
                }
                if (response.status == 200) {
                    return response.json()
                } else {
                    throw Error(response.statusText)
                }
          })
          .then(data => {
            this.setState({
                articles: data.data, DetailArticle: data.data[0], articlesInfo: data.columns, message: data.message
            });
          })
          .catch(error => {
            console.log(error);
                this.setState({ articles: [], articlesInfo: [] });
          })
    }

    createTask() {
        // Отправляем запрос на сервер для получения статей
        let query_url = ''
        try {
            query_url = variables.API_URL + this.createQuery()
        } catch(e) {
            alert('Не введен текст для поиска');
            return
        }
        fetch(query_url,
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${this.state.token}`,
                },
            })
            .then(response => response.json())
            .then(data => {
                this.setState({
                    full_query: data.task.full_query,
                    translation_stack: data.task.translation_stack,
                    short_query: data.task.query,
                    task: data.task,
                    count: data.task.count,
                    articles: [],
                    articlesInfo: []
                });
                alert("You query in queue? please wait to get result");
                this.getArticles()
            })
            .catch(error => {
                console.log(error);
                this.setState({ task: null });
            }
        )
    }

    getGraphData() {
        fetch(variables.API_URL + '/graphs',
            {
                headers: {
                    'Content-Type': 'application/json;charset=utf-8',
                    'Authorization': `Token ${variables.token}`,
                },
            }
          )
          .then(response => {
                console.log(response.status);
                if (response.status == 200) {
                    return response.json()
                } else {
                    throw Error(response.statusText)
                }
          })
          .then(data => {
            this.setState({
                graph: data.graph
            });
          })
          .catch(error => {
            console.log(error);
                this.setState({ graph: null});
          })
    }

    componentDidMount() {
        this.getArticles();
        console.log('start');
    }

    onSelectionChanged = () => {
        const selectedRows = this.gridRef.current.api.getSelectedRows();
        this.setState({DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null)})
    }

    changeQueryText = (e) => {
        this.setState({ queryText: e.target.value });
    }

    changeQueryStartDate = (e) => {
        this.setState({ queryStartDate: e.target.value });
    }

    changeQueryEndDate = (e) => {
        this.setState({ queryEndDate: e.target.value });
    }

    changeQueryGenders(gender) {
        if (this.state.queryGenders.has(gender)) {
            this.state.queryGenders.delete(gender)
        } else {
            this.state.queryGenders.add(gender)
        }
        this.setState({updateOr: !this.state.updateOr})

    }

    changeQueryTypes(type) {
        if (this.state.queryTypes.has(type)) {
            this.state.queryTypes.delete(type)
        } else {
            this.state.queryTypes.add(type)
        }
        this.setState({updateOr: !this.state.updateOr})
    }

    changeQueryOlds(old) {
        if (this.state.queryOlds.has(old)) {
            this.state.queryOlds.delete(old)
        } else {
            this.state.queryOlds.add(old)
        }
        this.setState({updateOr: !this.state.updateOr})
    }

    changeQueryOldsMany(olds) {
        if (!this.state.useAll) {
            for (let old of olds) {
                this.state.queryOlds.delete(old)
            }
        } else {
            for (let old of olds) {
                this.state.queryOlds.add(old)
            }
        }
        this.setState({useAll: !this.state.useAll})
    }

//    selectEvent() {
//        stabilized: () => {
//            if (network) { // Network will be set using getNetwork event from the Graph component
//                network.setOptions({ physics: false }); // Disable physics after stabilization
//                network.fit();
//            }
//        }
//    }

    render() {
        const {
            token,
            count,
            articles,
            DetailArticle,
            articlesInfo,
            translation_stack,
            full_query,
            short_query,
            message,

            queryText,
            queryEndDate,
            queryStartDate,

            graph,
        } = this.state;

        if (!token){
            return <Navigate push to="/login" />
        } else {
            return (
            <>
                <div className="container-fluid">
                    <div className="row">
                      <header className="navbar navbar-dark sticky-top bg-primary flex-md-nowrap p-0 shadow-sm">
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-2">
                            <button className="mt-1 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-controls="sidebar" aria-expanded="false" aria-label="Toggle navigation">
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                          <div className="col-md-10">
                            <a className="navbar-brand" href="#">SECHENOV AI-DATAMED</a>
                          </div>
                        </div>
                      </div>
                      <div className="col-md-8">
                        <input
                            className="form-control w-100"
                            id="search"
                            type="text"
                            name="search_field"
                            placeholder="Поисковый запрос"
                            value={queryText}
                            onChange={this.changeQueryText}
                            aria-label="Search" />
                      </div>
                      <div className="col-md-2">
                        <div className="row g-0">
                          <div className="col-md-10">
                            <ul className="navbar-nav px-3">
                              <li className="nav-item text-nowrap">
                                <input className="btn btn-primary" type="submit" value="Найти" onClick={() => this.createTask()}/>
                              </li>
                            </ul>
                          </div>
                          <div className="col-md-2">
                            <button className="mt-2 navbar-toggler collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#sidebar2" aria-controls="sidebar2" aria-expanded="false" aria-label="Toggle navigation">
                              <span className="icon-bar top-bar"></span>
                              <span className="icon-bar middle-bar"></span>
                              <span className="icon-bar bottom-bar"></span>
                            </button>
                          </div>
                        </div>
                      </div>
                      </header>
                    </div>
                </div>
                <main>
                    <div>
                        <div className="container-fluid">
                            <div className="row align-items-stretch b-height">

                            <aside id="sidebar" className="col-md-2 bg-light collapse show width mb-5 shadow-sm g-0">
                                <div className="accordion accordion-flush" id="accordionFlushExample">
                                  <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingOne">
                                      <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseOne" aria-expanded="false" aria-controls="flush-collapseOne">
                                        Дата публикации
                                      </button>
                                    </h2>
                                    <div id="flush-collapseOne" className="collapse show multi-collapse" aria-labelledby="flush-headingOne" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="mb-3">
                                          <label for="localdate">От : </label>
                                          <input
                                            type="date"
                                            id="d1"
                                            name="dateStart"
                                            value={queryStartDate}
                                            onChange={this.changeQueryStartDate}
                                          />
                                        </div>
                                        <div>
                                          <label for="localdate">До : </label>
                                          <input
                                            type="date"
                                            id="d2"
                                            name="dateStop"
                                            value={queryEndDate}
                                            onChange={this.changeQueryEndDate}
                                          />
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                  <div className="accordion-item">
                                    <h2 class="accordion-header" id="flush-headingThree">
                                      <button class="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseFour" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Тип статьи
                                      </button>
                                    </h2>
                                    <div id="flush-collapseFour" className="collapse show multi-collapse" aria-labelledby="flush-headingFour" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxClinicalTrial"
                                            name = "CheckBoxClinicalTrial"
                                            checked={this.state.queryTypes.has('clinical trial')}
                                            onChange={() => this.changeQueryTypes('clinical trial')}
                                          />
                                          <label className="form-check-label" for="CheckboxClinicalTrial">Clinical Trial</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMetaAnalysys"
                                            name = "CheckboxMetaAnalysys"
                                            checked={this.state.queryTypes.has('meta-analysis')}
                                            onChange={() => this.changeQueryTypes('meta-analysis')}
                                          />
                                          <label className="form-check-label" for="CheckboxMetaAnalysys">Meta Analysys</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxRandomizedControlledTrial"
                                            name ="CheckboxRandomizedControlledTrial"
                                            checked={this.state.queryTypes.has('randomized controlled trial')}
                                            onChange={() => this.changeQueryTypes('randomized controlled trial')}
                                          />
                                          <label className="form-check-label" for="CheckboxRandomizedControlledTrial">Randomized Controlled Trial</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxReview"
                                            name ="CheckboxReview"
                                            checked={this.state.queryTypes.has('review')}
                                            onChange={() => this.changeQueryTypes('review')}
                                          />
                                          <label className="form-check-label" for="CheckboxReview">Review</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxSystematicReview"
                                            name ="CheckboxSystematicReview"
                                            checked={this.state.queryTypes.has('systematic review')}
                                            onChange={() => this.changeQueryTypes('systematic review')}
                                          />
                                          <label className="form-check-label" for="CheckboxSystematicReview">Systematic Review</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxJournalArticle"
                                            name ="CheckboxJournalArticle"
                                            checked={this.state.queryTypes.has('journal article')}
                                            onChange={() => this.changeQueryTypes('journal article')}
                                          />
                                          <label className="form-check-label" for="CheckboxJournalArticle">Journal Article</label>
                                        </div>

                                      </div>
                                    </div>
                                  </div>
                                  <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingThree">
                                      <button className="accordion-button collapsed" data-target='#flush-collapseThree' type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseThree" aria-expanded="false" aria-controls="flush-collapseThree">
                                        Возраст пациента
                                      </button>
                                    </h2>
                                    <div id="flush-collapseThree" className="collapse show multi-collapse" aria-labelledby="flush-headingThree" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxChild18"
                                            checked={this.state.queryOlds.has('infant[mh]') & this.state.queryOlds.has('child[mh]') & this.state.queryOlds.has('adolescent[mh]')}
                                            onChange={() => this.changeQueryOldsMany(['infant[mh]','child[mh]', 'adolescent[mh]'])}
                                          />
                                          <label className="form-check-label" for="CheckboxChild18">Child: birth-18 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxNewborn"
                                            checked={this.state.queryOlds.has('infant, newborn[mh]')}
                                            onChange={() => this.changeQueryOlds('infant, newborn[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxNewborn">Newborn: birth-1 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxInfant023"
                                            checked={this.state.queryOlds.has('infant[mh]')}
                                            onChange={() => this.changeQueryOlds('infant[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxInfant023">Infant: birth-23 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxInfant123"
                                            checked={this.state.queryOlds.has('infant[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('infant[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxInfant123">Infant: 1-23 months</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxPreschool"
                                            checked={this.state.queryOlds.has('child, preschool[mh]')}
                                            onChange={() => this.changeQueryOlds('child, preschool[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxPreschool">Preschool Child: 2-5 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxChild612"
                                            checked={this.state.queryOlds.has('child[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('child[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxChild612">Child: 6-12 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdolescent"
                                            checked={this.state.queryOlds.has('adolescent[mh]')}
                                            onChange={() => this.changeQueryOlds('adolescent[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdolescent">Adolescent: 13-18 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdult19"
                                            checked={this.state.queryOlds.has('adult[mh]')}
                                            onChange={() => this.changeQueryOlds('adult[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdult19">Adult: 19+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxYoungAdult"
                                            checked={this.state.queryOlds.has('young adult[mh]')}
                                            onChange={() => this.changeQueryOlds('young adult[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxYoungAdult">Young Adult: 19-24 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAdult1944"
                                            checked={this.state.queryOlds.has('adult[mh:noexp]')}
                                            onChange={() => this.changeQueryOlds('adult[mh:noexp]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAdult1944">Adult: 19-44 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMiddle45"
                                            checked={this.state.queryOlds.has('middle aged[mh]') & this.state.queryOlds.has('aged[mh]')}
                                            onChange={() => this.changeQueryOldsMany(['middle aged[mh]', 'aged[mh]'])}
                                          />
                                          <label className="form-check-label" for="CheckboxMiddle45">Middle Aged: + Aged 45+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxMiddle4565"
                                            checked={this.state.queryOlds.has('middle aged[mh]')}
                                            onChange={() => this.changeQueryOlds('middle aged[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxMiddle4565">Middle Aged: 45-65 years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="CheckboxAged"
                                            checked={this.state.queryOlds.has('aged[mh]')}
                                            onChange={() => this.changeQueryOlds('aged[mh]')}
                                          />
                                          <label className="form-check-label" for="CheckboxAged">Aged: 65+ years</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            id="Checkbox80"
                                            checked={this.state.queryOlds.has('aged, 80 and over[mh]')}
                                            onChange={() => this.changeQueryOlds('aged, 80 and over[mh]')}
                                          />
                                          <label className="form-check-label" for="Checkbox80">80 and over: 80+ years</label>
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                    <div className="accordion-item">
                                    <h2 className="accordion-header" id="flush-headingTwo">
                                      <button className="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseTwo" aria-expanded="false" aria-controls="flush-collapseTwo">
                                        Пол
                                      </button>
                                    </h2>
                                    <div id="flush-collapseTwo" className="collapse show multi-collapse" aria-labelledby="flush-headingTwo" data-bs-target="#accordionFlushExample">
                                      <div className="accordion-body">
                                        <div className="form-check form-check-inline">
                                          <input
                                            type="checkbox"
                                            id="CheckboxMan"
                                            name = "CheckboxMan"
                                            className="form-check-input"
                                            checked={this.state.queryGenders.has('man')}
                                            onChange={() => this.changeQueryGenders('man')}
                                          />
                                          <label className="form-check-label" for="CheckboxMan">Мужчина</label>
                                        </div>
                                        <div className="form-check form-check-inline">
                                          <input
                                            className="form-check-input"
                                            type="checkbox"
                                            name = "CheckboxWomen"
                                            id="CheckboxWomen"
                                            checked={this.state.queryGenders.has('woman')}
                                            onChange={() => this.changeQueryGenders('woman')}
                                          />
                                          <label className="form-check-label" for="CheckboxWomen">Женщина</label>
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                </div>
                            </aside>

                            <section class="col shadow p-4" style={{backgroundColor: "#fff"}}>
                              <div class="accordion accordion-flush" id="accordion">
                                <div class="accordion-item">
                                  <h2 class="accordion-header" id="">
                                    <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#flush-collapseSeven" aria-expanded="false" aria-controls="flush-collapseSeven">
                                      По запросу найдено { count } источников.
                                    </button>
                                  </h2>
                                  <div id="flush-collapseSeven" class="collapse multi-collapse" aria-labelledby="flush-headingSeven" data-bs-target="#accordionFlushExample">
                                    <div class="accordion-body">
                                      <p class="pb-2 mb-3 border-bottom"> Запрос {short_query} .</p>
                                      <p class="pb-2 mb-3 border-bottom"> Запрос автоматически расширен до следующего вида - {full_query}.</p>
                                      <p class="pb-2 mb-3 border-bottom"> Служебная информация для анализа : {translation_stack}.</p>
                                    </div>
                                  </div>
                                </div>
                              </div>
                              <div>
                              <div class="bd-example">
                                <ul class="nav nav-pills mb-3" id="myTab" role="tablist">
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link active" id="home-tab" data-bs-toggle="tab" data-bs-target="#home" type="button" role="tab" aria-controls="home" aria-selected="false">Тематическое описание коллекции</button>
                                  </li>
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link" id="profile-tab" data-bs-toggle="tab" data-bs-target="#profile" type="button" role="tab" aria-controls="profile" aria-selected="true">Анализ</button>
                                  </li>
                                  <li class="nav-item" role="presentation">
                                    <button class="nav-link" id="contact-tab" data-bs-toggle="tab" data-bs-target="#contact" type="button" role="tab" aria-controls="contact" aria-selected="false" onClick={() => this.getGraphData()}>Схема</button>
                                  </li>
                                </ul>
                                <div class="tab-content" id="myTabContent">
                                  <div class="tab-pane fade active show" id="home" role="tabpanel" aria-labelledby="home-tab">
                                    <div class="container-fluid g-0">
                                        <div className="ag-theme-alpine" style={{height: 700}}>
                                            <AgGridReact
                                                ref={this.gridRef}
                                                rowData={articles}
                                                columnDefs={articlesInfo}
                                                pagination={true}
                                                rowSelection={'single'}
                                                onSelectionChanged={this.onSelectionChanged}
                                            >
                                            </AgGridReact>
                                        </div>
                                    </div>
                                  </div>
                                  <div class="tab-pane fade" id="profile" role="tabpanel" aria-labelledby="profile-tab">
                                    <div class="container-fluid g-0" style={{height: "700px"}}>
                                      <div class="row">
                                        <div class="col-md-4">
                                          <div class="col">
                                            <div class="card p-4" style={{height: "475px", overflow: "scroll"}}>
                                              <ul class="treeline" id="list">
                                                <li>Список 1
                                                  <ul>
                                                    <li>Подраздел 1
                                                      <ul>
                                                        <li>Пункт 1</li>
                                                        <li>Пункт 1</li>
                                                        <li>Пункт 1</li>
                                                        <li>Пункт 1</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 2
                                                      <ul>
                                                        <li>Пункт 2</li>
                                                        <li>Пункт 2</li>
                                                        <li>Пункт 3</li>
                                                      </ul>
                                                    </li>
                                                  </ul>
                                                </li>
                                                <li>Список 2
                                                  <ul>
                                                    <li>Подраздел 3
                                                      <ul>
                                                        <li>Пункт 3</li>
                                                        <li>Пункт 3</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 4
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 5
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 6
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 7
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 8
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                    <li>Подраздел 9
                                                      <ul>
                                                        <li>Пункт 4</li>
                                                        <li>Пункт 4</li>
                                                      </ul>
                                                    </li>
                                                  </ul>
                                                </li>
                                              </ul>
                                            </div>
                                          </div>
                                          <div class="col my-4">
                                            <div class="card" style={{height: "200px", overflow: "scroll"}}>
                                              <div class="tag-container p-4">
                                                <div class="tag-cloud">
                                                  <div class="cloud-item" tag-value="3">Aenean</div>
                                                  <div class="cloud-item" tag-value="8">Vulputate</div>
                                                  <div class="cloud-item" tag-value="4">Eleifend</div>
                                                  <div class="cloud-item" tag-value="1">Harum</div>
                                                  <div class="cloud-item" tag-value="9">Dapibus</div>
                                                  <div class="cloud-item" tag-value="6">Rerum</div>
                                                  <div class="cloud-item" tag-value="4">Deleniti</div>
                                                  <div class="cloud-item" tag-value="5">Cupiditate</div>
                                                  <div class="cloud-item" tag-value="1">Nisi</div>
                                                  <div class="cloud-item" tag-value="8">Ratione</div>
                                                  <div class="cloud-item" tag-value="8">Incidunt</div>
                                                  <div class="cloud-item" tag-value="2">Odio</div>
                                                </div>
                                              </div>
                                            </div>
                                          </div>
                                        </div>
                                        <div class="col-md-8">
                                          <div class="card p-4" style={{height: "700px", overflow: "scroll"}}>
                                          <p>Cras dapibus. Vivamus elementum semper nisi. Aenean vulputate eleifend tellus. Aenean leo ligula, porttitor eu, consequat vitae, eleifend ac, enim. Aliquam lorem ante, dapibus in, viverra quis, feugiat a, tellus. Phasellus viverra nulla ut metus varius laoreet. Quisque rutrum. Aenean imperdiet. Etiam ultricies nisi vel augue. Curabitur ullamcorper ultricies nisi. Nam eget dui. Etiam rhoncus. Maecenas tempus, tellus eget condimentum rhoncus, sem quam semper libero, sit amet adipiscing sem neque sed ipsum. Nam quam nunc, blandit vel, luctus pulvinar, hendrerit id, lorem. Maecenas nec odio et ante tincidunt tempus. Donec vitae sapien ut libero venenatis faucibus. Nullam quis ante.
                                          Cras dapibus. Vivamus elementum semper nisi. Aenean vulputate eleifend tellus. Aenean leo ligula, porttitor eu, consequat vitae, eleifend ac, enim. Aliquam lorem ante, dapibus in, viverra quis, feugiat a, tellus. Phasellus viverra nulla ut metus varius laoreet. Quisque rutrum. Aenean imperdiet. Etiam ultricies nisi vel augue. Curabitur ullamcorper ultricies nisi. Nam eget dui. Etiam rhoncus. Maecenas tempus, tellus eget condimentum rhoncus, sem quam semper libero, sit amet adipiscing sem neque sed ipsum. Nam quam nunc, blandit vel, luctus pulvinar, hendrerit id, lorem. Maecenas nec odio et ante tincidunt tempus. Donec vitae sapien ut libero venenatis faucibus. Nullam quis ante.
                                          Cras dapibus. Vivamus elementum semper nisi. Aenean vulputate eleifend tellus. Aenean leo ligula, porttitor eu, consequat vitae, eleifend ac, enim. Aliquam lorem ante, dapibus in, viverra quis, feugiat a, tellus. Phasellus viverra nulla ut metus varius laoreet. Quisque rutrum. Aenean imperdiet. Etiam ultricies nisi vel augue. Curabitur ullamcorper ultricies nisi. Nam eget dui. Etiam rhoncus. Maecenas tempus, tellus eget condimentum rhoncus, sem quam semper libero, sit amet adipiscing sem neque sed ipsum. Nam quam nunc, blandit vel, luctus pulvinar, hendrerit id, lorem. Maecenas nec odio et ante tincidunt tempus. Donec vitae sapien ut libero venenatis faucibus. Nullam quis ante.
                                          Cras dapibus. Vivamus elementum semper nisi. Aenean vulputate eleifend tellus. Aenean leo ligula, porttitor eu, consequat vitae, eleifend ac, enim. Aliquam lorem ante, dapibus in, viverra quis, feugiat a, tellus. Phasellus viverra nulla ut metus varius laoreet. Quisque rutrum. Aenean imperdiet. Etiam ultricies nisi vel augue. Curabitur ullamcorper ultricies nisi. Nam eget dui. Etiam rhoncus. Maecenas tempus, tellus eget condimentum rhoncus, sem quam semper libero, sit amet adipiscing sem neque sed ipsum. Nam quam nunc, blandit vel, luctus pulvinar, hendrerit id, lorem. Maecenas nec odio et ante tincidunt tempus. Donec vitae sapien ut libero venenatis faucibus. Nullam quis ante.</p>
                                          </div>
                                        </div>
                                      </div>
                                    </div>
                                  </div>
                                  <div class="tab-pane fade" id="contact" role="tabpanel" aria-labelledby="contact-tab">
                                    <div class="container-fluid g-0">
                                      <div id="mynetwork" style={{width: "100%", height: "600px"}}>
                                      {graph?
                                        <Graph
                                          graph={graph}
                                          zoomView={true}
//                                          events={this.selectEvent}
                                          key={uuidv4()}
                                          options={{
                                            layout: {
                                                hierarchical: false,
                                            },
                                            edges: {
                                              color: "#000000"
                                            }
                                          }}
//                                          getNetwork={network => {
//                                            setNetwork(network);
//                                          }}
                                        />
                                        :null}
                                      </div>
                                    </div>
                                  </div>
                                </div>
                              </div>
                              </div>
                            </section>

                            <aside id="sidebar2" class="col-md-4 bg-light collapse show width mb-5 shadow">
                              <h1 class="h2 pt-3 pb-2 mb-3 border-bottom">Подробности</h1>
                              <nav class="small" id="toc">
                                {DetailArticle?
                                    <div class="card mb-3">
                                        <div class="card-body">
                                          <a href= { DetailArticle.url } class="card-title link-primary text-decoration-none h5"> { DetailArticle.titl } </a>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Авторы :  { DetailArticle.auth } </p>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text">Аннотация :  </p>
                                          <p class="card-text"> { DetailArticle.tiab } </p>
                                          <p class="card-text">---------------------------------- </p>
                                          <p class="card-text"><small class="text-success">Дата публикации : { DetailArticle.pdat } </small></p>
                                          <p class="card-text"><small class="text-success">Издание : { DetailArticle.jour }</small></p>
                                          <p class="card-text"><small class="text-success">Вид публикации : { DetailArticle.pt }</small></p>
                                          <p class="card-text"><small class="text-success">Страна : { DetailArticle.pl } </small></p>
                                          <p class="card-text"><small class="text-success">{ DetailArticle.mesh } </small></p>
                                        </div>
                                      </div>
                                :null}
                              </nav>
                            </aside>

                            </div>
                        </div>
                    </div>
                </main>
            </>
            )
        }
    }
}